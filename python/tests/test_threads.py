"""Threads: the GIL is released, separate folders run in parallel, and the engine's
internal rayon parallelism gives the same results at any thread count.

The engine parallelises its inner loops with rayon (one global pool per process,
sized by ``RAYON_NUM_THREADS``). It keeps the box geometry in a thread-local
variable that every parallel task re-syncs; the tests with the triclinic box of
``EPSR-test-ex14`` check this through the wrapper, including two simulations
with different boxes sharing the pool at the same time.
"""

from __future__ import annotations

import os
import subprocess
import sys
import threading
import time
from collections.abc import Callable
from pathlib import Path

import pytest

import epsr

from epsr_cases import (
    CASE_FILES_TO_SKIP,
    TEMPLATE,
    TEMPLATE_INPUT_FILE,
    TEST_PROBLEMS,
    compare_dirs,
    copy_files,
)

#: The EPSR-test system in a triclinic (non-cubic) box.
TRICLINIC_CASE = TEST_PROBLEMS / "200DMSO400Water" / "EPSR-test-ex14"
BOX_SHAPES = {"cubic": TEMPLATE, "triclinic": TRICLINIC_CASE}

#: Runs one EPSR invocation in a fresh process, whose rayon pool size is set by the environment.
RUN_ONE_SIMULATION = "import sys, epsr; epsr.Simulation(sys.argv[1], sys.argv[2]).run()"

#: Iterations of a pure-Python loop that must complete while the engine runs in another thread.
MINIMUM_MAIN_THREAD_ITERATIONS = 1000


def _run_in_threads(work: Callable[[int], None], n_threads: int) -> None:
    threads = [
        threading.Thread(target=work, args=(index,)) for index in range(n_threads)
    ]
    for thread in threads:
        thread.start()
    for thread in threads:
        thread.join()


def test_run_releases_the_gil(project: epsr.Project) -> None:
    simulation = project.simulation()
    run_finished = threading.Event()

    def run_simulation() -> None:
        simulation.run()
        run_finished.set()

    engine_thread = threading.Thread(target=run_simulation)
    main_thread_iterations = 0
    engine_thread.start()
    while not run_finished.is_set():  # pure-Python work on the main thread
        main_thread_iterations += 1
        time.sleep(0)
    engine_thread.join()

    assert simulation.completed
    # If the engine held the GIL for the whole run, this loop could not advance.
    assert main_thread_iterations > MINIMUM_MAIN_THREAD_ITERATIONS


def test_parallel_runs_match_sequential(template: Path, tmp_path: Path) -> None:
    sequential_folder = copy_files(template, tmp_path / "sequential")
    sequential_summary = epsr.Simulation(sequential_folder).run()

    n_parallel = 3
    parallel_folders = [
        copy_files(template, tmp_path / f"parallel{index}")
        for index in range(n_parallel)
    ]
    parallel_summaries: dict[int, epsr.Summary] = {}

    def run_one(index: int) -> None:
        parallel_summaries[index] = epsr.Simulation(parallel_folders[index]).run()

    _run_in_threads(run_one, n_parallel)

    assert all(
        parallel_summaries[index] == sequential_summary for index in range(n_parallel)
    )
    sequential_pcof = (sequential_folder / "200DMSO400Water.pcof").read_bytes()
    for folder in parallel_folders:
        assert (folder / "200DMSO400Water.pcof").read_bytes() == sequential_pcof


def test_refine_many(template: Path, tmp_path: Path) -> None:
    projects = [
        epsr.Project.from_template(template, tmp_path / f"project{index}")
        for index in range(2)
    ]
    histories = epsr.refine_many(projects, 2)
    assert len(histories) == 2
    assert all(len(history) == 2 for history in histories)
    assert list(histories[0].r_factor) == list(histories[1].r_factor)
    with pytest.raises(ValueError, match="more than once"):
        epsr.refine_many([projects[0], projects[0].folder], 1)


def test_tools_from_threads(template: Path, tmp_path: Path) -> None:
    """Tool runs from several threads are serialised safely and give the shipped weights."""
    n_threads = 4
    folders = [
        copy_files(template, tmp_path / f"tool{index}") for index in range(n_threads)
    ]
    errors: list[BaseException] = []

    def run_nwts(index: int) -> None:
        try:
            epsr.tools.nwts(folders[index], "DMSODWaterH.NWTS.dat")
        except BaseException as error:  # reported by the assertion below
            errors.append(error)

    _run_in_threads(run_nwts, n_threads)

    assert not errors
    expected_weights = (template / "DMSODWaterH.NWTStot.wts").read_bytes()
    for folder in folders:
        assert (folder / "DMSODWaterH.NWTStot.wts").read_bytes() == expected_weights
    assert Path.cwd() not in folders  # the working directory was restored


def _copy_case(source: Path, destination: Path) -> Path:
    if not source.exists():
        pytest.skip(f"{source} not found")
    return copy_files(source, destination, CASE_FILES_TO_SKIP)


def _run_with_rayon_threads(folder: Path, rayon_threads: int) -> None:
    """One wrapper run in a new process: the rayon pool size is fixed per process."""
    environment = {**os.environ, "RAYON_NUM_THREADS": str(rayon_threads)}
    subprocess.run(
        [sys.executable, "-c", RUN_ONE_SIMULATION, str(folder), TEMPLATE_INPUT_FILE],
        env=environment,
        check=True,
        capture_output=True,
    )


@pytest.mark.parametrize("box_shape", list(BOX_SHAPES))
def test_results_do_not_depend_on_rayon_thread_count(
    box_shape: str, tmp_path: Path
) -> None:
    """The wrapper writes byte-identical files with 1 and with 4 engine threads."""
    source = BOX_SHAPES[box_shape]
    single_thread_folder = _copy_case(source, tmp_path / "one_thread")
    multi_thread_folder = _copy_case(source, tmp_path / "four_threads")

    _run_with_rayon_threads(single_thread_folder, 1)
    _run_with_rayon_threads(multi_thread_folder, 4)

    comparison = compare_dirs(single_thread_folder, multi_thread_folder)
    assert comparison.all_identical, comparison.report()


def test_concurrent_cubic_and_triclinic_runs_match_sequential(tmp_path: Path) -> None:
    """Two simulations with different boxes share the rayon pool without mixing geometries."""
    sequential_folders = {
        shape: _copy_case(source, tmp_path / f"sequential_{shape}")
        for shape, source in BOX_SHAPES.items()
    }
    concurrent_folders = {
        shape: _copy_case(source, tmp_path / f"concurrent_{shape}")
        for shape, source in BOX_SHAPES.items()
    }
    for folder in sequential_folders.values():
        epsr.Simulation(folder, TEMPLATE_INPUT_FILE).run()

    shapes = list(BOX_SHAPES)
    errors: list[BaseException] = []

    def run_one(index: int) -> None:
        try:
            for _ in range(2):  # two chained runs keep both simulations busy together
                epsr.Simulation(
                    concurrent_folders[shapes[index]], TEMPLATE_INPUT_FILE
                ).run()
        except BaseException as error:  # reported by the assertion below
            errors.append(error)

    for folder in sequential_folders.values():  # the second sequential run
        epsr.Simulation(folder, TEMPLATE_INPUT_FILE).run()
    _run_in_threads(run_one, len(shapes))

    assert not errors
    for shape in shapes:
        comparison = compare_dirs(sequential_folders[shape], concurrent_folders[shape])
        assert comparison.all_identical, f"{shape} box:\n{comparison.report()}"
