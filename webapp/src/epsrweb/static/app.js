/* Results page: tabs, the live event stream, Overleaf publishing and chat. */
(function () {
  "use strict";

  var root = document.getElementById("runroot");
  if (!root) return;

  var runId = root.dataset.runId;
  var csrf = root.dataset.csrf;
  var TAB_KEY = "epsrweb.tab." + runId;

  /* ------------------------------------------------------------------ tabs */

  var tabs = Array.prototype.slice.call(document.querySelectorAll(".tab"));

  function show(panelId) {
    tabs.forEach(function (tab) {
      var on = tab.dataset.panel === panelId;
      tab.classList.toggle("active", on);
      tab.setAttribute("aria-selected", on ? "true" : "false");
    });
    document.querySelectorAll(".panel").forEach(function (panel) {
      panel.classList.toggle("active", panel.id === panelId);
    });
    try { sessionStorage.setItem(TAB_KEY, panelId); } catch (e) { /* private mode */ }
  }

  tabs.forEach(function (tab) {
    tab.addEventListener("click", function () { show(tab.dataset.panel); });
  });

  try {
    var remembered = sessionStorage.getItem(TAB_KEY);
    if (remembered && document.getElementById(remembered)) show(remembered);
  } catch (e) { /* ignore */ }

  /* ------------------------------------------------------- the event stream */

  var log = document.getElementById("log");
  var pill = document.getElementById("statuspill");
  var status = root.dataset.status;
  var live = status === "queued" || status === "running";

  function append(line) {
    if (!log || !line) return;
    var atBottom = log.scrollHeight - log.scrollTop - log.clientHeight < 40;
    log.textContent += line + "\n";
    if (atBottom) log.scrollTop = log.scrollHeight;
  }

  if (live) {
    if (log) log.scrollTop = log.scrollHeight;
    var source = new EventSource("/runs/" + encodeURIComponent(runId) + "/events");

    source.onmessage = function (message) {
      try {
        append(JSON.parse(message.data).line);
      } catch (e) { /* a malformed frame is not worth a broken page */ }
    };

    source.addEventListener("end", function (message) {
      source.close();
      var finalStatus = "done";
      try { finalStatus = JSON.parse(message.data).status || "done"; } catch (e) { /* keep */ }
      if (pill) {
        pill.textContent = finalStatus;
        pill.className = "pill " + finalStatus;
      }
      append("-- " + finalStatus + ". Reloading to show the results.");
      /* The figures, the analysis and the downloads are rendered server-side,
         so the page reloads once the job is over; the open tab is restored. */
      window.setTimeout(function () { window.location.reload(); }, 900);
    });

    source.onerror = function () {
      /* Let EventSource retry on its own; speak up only if it gave up. */
      if (source.readyState === EventSource.CLOSED) {
        append("-- the progress stream closed; reload the page to catch up.");
      }
    };
  }

  /* -------------------------------------------------------- Overleaf push */

  var form = document.getElementById("figform");
  var pushButton = document.getElementById("pushbtn");
  var pushResult = document.getElementById("pushresult");

  function selectedFigures() {
    if (!form) return [];
    return Array.prototype.slice
      .call(form.querySelectorAll('input[name="figures"]:checked'))
      .map(function (box) { return box.value; });
  }

  function selectedFormat() {
    var picked = form && form.querySelector('input[name="figure_format"]:checked');
    return picked ? picked.value : "pdf";
  }

  function escapeHtml(text) {
    var holder = document.createElement("div");
    holder.textContent = text;
    return holder.innerHTML;
  }

  function report(kind, html) {
    if (!pushResult) return;
    pushResult.hidden = false;
    pushResult.className = "pushresult strip " + kind;
    pushResult.innerHTML = html;
  }

  function asJson(response) {
    return response.json().then(function (data) {
      return { ok: response.ok, data: data || {} };
    });
  }

  if (pushButton) {
    pushButton.addEventListener("click", function () {
      var figures = selectedFigures();
      if (!figures.length) {
        report("warn", "Tick at least one figure first.");
        return;
      }
      var original = pushButton.textContent;
      pushButton.disabled = true;
      pushButton.textContent = "Pushing...";
      report("", "Cloning the project, copying " + figures.length + " figure(s), pushing...");

      fetch("/runs/" + encodeURIComponent(runId) + "/overleaf", {
        method: "POST",
        headers: { "Content-Type": "application/json", "X-CSRF-Token": csrf },
        body: JSON.stringify({ figures: figures, format: selectedFormat() })
      })
        .then(asJson)
        .then(function (result) {
          var data = result.data;
          if (!result.ok || data.ok === false) {
            report("bad", "<strong>Not pushed.</strong> " +
              escapeHtml(data.detail || "Unknown error."));
            return;
          }
          var html = "<strong>Pushed.</strong> " + escapeHtml(data.detail || "");
          if (data.commit) {
            html += " Commit <code>" + escapeHtml(data.commit) + "</code>.";
          }
          if (data.snippets && data.snippets.length) {
            html += "<p>Paste into your document:</p><pre class='code'>" +
              escapeHtml(data.snippets.join("\n")) + "</pre>";
          }
          report("ok", html);
        })
        .catch(function (error) {
          report("bad", "<strong>Not pushed.</strong> " + escapeHtml(String(error)));
        })
        .then(function () {
          pushButton.disabled = false;
          pushButton.textContent = original;
        });
    });
  }

  if (form) {
    form.addEventListener("submit", function (event) {
      if (!selectedFigures().length) {
        event.preventDefault();
        report("warn", "Tick at least one figure first.");
      }
    });
  }

  /* ------------------------------------------------------------------ chat */

  var chatForm = document.getElementById("chatform");
  var chatInput = document.getElementById("chatinput");
  var chatLog = document.getElementById("chatlog");

  function bubble(text, cls) {
    var div = document.createElement("div");
    div.className = "bubble " + cls;
    div.textContent = text;
    chatLog.appendChild(div);
    div.scrollIntoView({ block: "nearest" });
    return div;
  }

  if (chatForm) {
    chatForm.addEventListener("submit", function (event) {
      event.preventDefault();
      var message = (chatInput.value || "").trim();
      if (!message) return;
      chatInput.value = "";
      chatInput.disabled = true;
      bubble(message, "you");
      var waiting = bubble("thinking, and querying the run's results...", "pending");

      fetch("/runs/" + encodeURIComponent(runId) + "/chat", {
        method: "POST",
        headers: { "Content-Type": "application/json", "X-CSRF-Token": csrf },
        body: JSON.stringify({ message: message })
      })
        .then(asJson)
        .then(function (result) {
          waiting.remove();
          bubble(result.ok ? (result.data.answer || "(no answer)")
                           : (result.data.detail || "The agent could not answer."),
                 result.ok ? "them" : "pending");
        })
        .catch(function (error) {
          waiting.remove();
          bubble(String(error), "pending");
        })
        .then(function () {
          chatInput.disabled = false;
          chatInput.focus();
        });
    });
  }
})();
