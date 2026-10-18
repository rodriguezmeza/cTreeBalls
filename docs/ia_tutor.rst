Interactive scientific tutor
============================

Start the tutor from a built source checkout::

   ./cballs options=IA-tutor

If ``cballs`` is on your executable path, omit ``./``. This separate mode
starts Python before catalogue allocation or MPI initialization. It does not
start a correlation calculation. Python 3.9 or newer is required; the default
mode uses only its standard library and needs no model, network, API key or
Cython extension. Direct Python access is also available::

   python3 scripts/ia_tutor.py
   ./cballs options=IA-tutor --ask "Does bin theta override tree theta?"
   ./cballs options=IA-tutor --ask "What is zero tolerance?" --json

The default is an **offline documentation tutor**, not a generative language
model. It searches current documentation and supplies reviewed scientific
explanations for the questions covered by its local knowledge. Unknown
questions receive relevant excerpts or an explicit statement that evidence is
missing. Answers include file references. An optional local model can compose
more flexible explanations from retrieved passages.

A conversation
--------------

For example::

   You> /context field=shear geometry=sphere statistic=3pcf threads=16
   You> Which engine should I compare for DESY3?
   You> Does CBALLS_SHEAR_BIN_THETA overwrite tree-theta?
   You> What is the role of CBALLS_SHEAR_PROFILE?
   You> Is it the same as verbose?
   You> /sources
   You> /topic measured-penalty
   You> /quit

The session retains field, geometry, statistic and recent question context.
Offline follow-up handling is deliberately limited; naming the parameter or
using ``/topic`` resolves ambiguous questions. ``/reset`` clears this context.
Approximation settings are never presented as universal relative-error bounds.
Historical timings and catalogue comparisons are labelled as examples; an
unmeasured penalty is not zero error.

Useful commands
---------------

* ``/help`` lists all commands; ``/topics`` lists reviewed explanations.
* ``/topic defaults`` explains defaults and explicit environment overrides.
* ``/search WORDS`` searches source passages; ``/source N`` shows the cited
  lines from the last answer.
* ``/engines shear`` checks available methods against this executable's
  ``print-search-methods`` output; ``/build`` reads its ``make-info`` output.
  A missing executable is reported as unverified availability.
* ``/catalog "PATH"`` inspects FITS or NPZ headers, suggests candidate fields,
  and updates catalogue context. It does not load the survey arrays or validate
  coordinate units, weights, masks or data quality. FITS inspection needs the
  optional ``astropy`` package; NPZ header inspection uses the standard library.
* ``/recommend`` compares methods compatible with the selected field and
  geometry. Specify ``statistic=2pcf``, ``3pcf`` or ``both`` as appropriate.
* ``/run "PATH"`` reads saved JSON run information or an
  ``approximation-penalty.json`` and explains its measured error table.
* ``/reindex`` refreshes the local document database.
* ``/save "PATH"`` saves the conversation as JSON, refusing to overwrite an
  existing file. Questions are not saved automatically. ``--transcript PATH``
  explicitly saves on normal exit; an interrupted session is not saved.
* ``/quit``, Ctrl-D or Ctrl-C exits.

Catalogue and run inspection are explicit local reads. The tutor does not
execute suggested analysis commands. Its only native subprocesses are fixed
``make-info`` and ``print-search-methods`` information queries.

Knowledge and updates
---------------------

The index includes the main README, documentation, public benchmark/regression
README files, enabled-addon README files, the engine capability table, and
``tutor/notes``. Reviewed explanations in ``tutor/knowledge/questions.json``
cover shear, kappa, Lyman-alpha, timing, approximation penalties, count multipoles,
normalization, binning, MPI, plotting and engine selection. They incorporate
scientific questions visible during development of this feature. The tutor
does not automatically access other chats, private catalogues, or directories
outside the selected checkout.

A SQLite cache is built at ``build/ia-tutor/knowledge.sqlite3``. Content hashes
refresh it when indexed source files change, are added or removed. An
unwritable or unusable cache falls back to an in-memory index. Use
``--database :memory:`` to avoid a persistent cache. ``--reindex`` forces a
rebuild; it does not certify old scientific explanations against changed
sources. Saved answers carry hashes of their reviewed references, and changed
references cause a visible warning with current documentation excerpts.

Add new scientific Markdown or reStructuredText notes under ``tutor/notes``.
They will be indexed at the next startup or ``/reindex``. See
``tutor/README.md`` for the reviewed question format. Keep numerical examples
explicit about estimator, input sample, approximation controls, error metric
and timing scope.

Optional local AI
-----------------

With a local Ollama service already running and a model already installed::

   ./cballs options=IA-tutor --backend ollama --model YOUR_INSTALLED_MODEL

The interface follows the official
`Ollama chat API <https://docs.ollama.com/api/chat>`_. Model installation and
service startup are separate user actions. This feature never downloads a
model. The default endpoint is ``http://127.0.0.1:11434``; ``--model-url`` accepts
loopback addresses only, including a local SSH tunnel. HTTP proxy environment
settings and redirects are not used. Questions, session context and selected
source excerpts are sent to that explicitly selected service. Treat the local
service and chosen model as part of your own software environment.

Model answers are labelled as local-model synthesis. Citations identify the
supplied evidence; they do not automatically verify a generated claim. If the
service is unavailable or its answer cannot be read, the tutor displays its
offline answer instead. ``--timeout`` controls the request timeout. Use
``--backend offline`` to explicitly disable model calls.

Configuration and deployment
----------------------------

The native launcher locates tutor resources near its executable, then at the
source directory recorded during compilation. A copied binary alone is not a
complete tutor installation. For a relocated checkout or a different Python::

   CBALLS_TUTOR_ROOT=/path/to/cTreeBalls \
   CBALLS_TUTOR_PYTHON=/path/to/python3 cballs options=IA-tutor

``CBALLS_TUTOR_PYTHON`` is one executable path, without shell arguments.
``CBALLS_TUTOR_ROOT`` selects the source checkout. Python options include
``--root`` and ``--cballs`` for explicit resource and build inspection paths.
``CBALLS_TUTOR_DATABASE``, ``CBALLS_TUTOR_BACKEND``, ``CBALLS_TUTOR_MODEL`` and
``CBALLS_TUTOR_OLLAMA_URL`` provide persistent preferences; command-line values
win. ``--context "field=kappa geometry=sphere statistic=3pcf"`` supplies initial
context. ``--help`` lists all arguments.

Use one terminal process for a conversation, including on a cluster. Do not
launch the interactive tutor across MPI ranks or Slurm tasks. The launcher
rejects detected multi-process launches. Combine neither ``IA-tutor`` with
other native ``options`` nor calculation parameters with tutor arguments.

Validation
----------

Run ``make test-ia-tutor`` in the selected build profile. Tests cover native
launching and argument forwarding, offline retrieval and follow-ups, document
refresh and stale-answer warnings, catalogue headers, saved penalties,
transcripts, and the optional model protocol with a simulated service.
