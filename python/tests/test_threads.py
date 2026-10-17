"""The GIL is released while the engine runs; separate folders can run in parallel."""

from __future__ import annotations

import threading
import time
from collections.abc import Callable
from pathlib import Path

import pytest

import epsr

from epsr_cases import copy_files

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
