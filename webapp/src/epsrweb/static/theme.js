/* Light / dark mode.

   Loaded in <head> without `defer`, so a stored choice is applied before the page
   is painted (no flash of the wrong theme). Without a stored choice the page
   follows the system setting through CSS alone; this script only adds the
   header button, which is hidden when scripts do not run. */
(function () {
  "use strict";

  var KEY = "epsrweb.theme";
  var root = document.documentElement;

  try {
    var stored = window.localStorage.getItem(KEY);
    if (stored === "light" || stored === "dark") root.setAttribute("data-theme", stored);
  } catch (e) { /* storage blocked: follow the system */ }

  var media = window.matchMedia ? window.matchMedia("(prefers-color-scheme: dark)") : null;

  function current() {
    var chosen = root.getAttribute("data-theme");
    if (chosen === "light" || chosen === "dark") return chosen;
    return media && media.matches ? "dark" : "light";
  }

  function render(button) {
    var next = current() === "dark" ? "light" : "dark";
    button.querySelector("[data-label]").textContent =
      next === "dark" ? "Dark mode" : "Light mode";
    button.setAttribute("aria-label", "Switch to " + next + " mode");
  }

  document.addEventListener("DOMContentLoaded", function () {
    var button = document.getElementById("theme-toggle");
    if (!button) return;
    button.hidden = false;
    render(button);
    button.addEventListener("click", function () {
      var next = current() === "dark" ? "light" : "dark";
      root.setAttribute("data-theme", next);
      try { window.localStorage.setItem(KEY, next); } catch (e) { /* not remembered */ }
      render(button);
    });
    if (media) {
      var onChange = function () { render(button); };
      if (media.addEventListener) media.addEventListener("change", onChange);
      else if (media.addListener) media.addListener(onChange);
    }
  });
})();
