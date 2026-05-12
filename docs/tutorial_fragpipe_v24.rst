FragPipe v24 DDA and DIA Tutorial
=================================

This tutorial shows how to prepare FragPipe v24.0 output for FLiPPR, then
run FLiPPR on DDA or DIA label-free data. It assumes that FragPipe has already
searched and quantified the raw mass spectrometry files. FLiPPR starts from the
FragPipe output directory.

FLiPPR is built for LiP-MS experiments with an optional matched TrP
normalization experiment. In this guide:

- ``LiP`` means the limited-proteolysis experiment.
- ``TrP`` means the trypsin-only experiment used for protein-level
  normalization.
- A ``process`` is one control-vs-test comparison inside a study.

The examples use two conditions, ``CTRL`` and ``DRUG``, each with three
replicates.

Inputs FLiPPR expects
---------------------

FLiPPR reads different FragPipe output files depending on the acquisition mode.

For DDA, the FragPipe output directory must contain:

.. code-block:: text

    combined_ion.tsv
    combined_protein.tsv
    experiment_annotation.tsv

For DIA, the FragPipe output directory must contain:

.. code-block:: text

    ion.tsv
    dia-quant-output/report.pr_matrix.tsv
    dia-quant-output/report.pg_matrix.tsv
    experiment_annotation.tsv

The DDA files come from FragPipe's Philosopher and IonQuant reports. The DIA
matrix files come from the DIA-NN quantification step run by FragPipe. See the
`FragPipe output guide`_ for the full list of reports and column meanings.

Sample naming
-------------

FLiPPR expects replicate intensity columns that follow the sample names in
``experiment_annotation.tsv``. Use stable sample names with the replicate number
at the end:

.. code-block:: text

    CTRL_1
    CTRL_2
    CTRL_3
    DRUG_1
    DRUG_2
    DRUG_3

When calling ``Study.add_process()``, pass the condition prefix without the
replicate suffix:

.. code-block:: python

    study.add_process(
        pid="drug",
        lip_ctrl="CTRL",
        lip_test="DRUG",
        n_rep=3,
    )

This tells FLiPPR to look for ``CTRL_1 Intensity``, ``CTRL_2 Intensity``,
``CTRL_3 Intensity``, ``DRUG_1 Intensity``, and so on.

For uneven replicate counts, pass a tuple:

.. code-block:: python

    study.add_process("drug", "CTRL", "DRUG", n_rep=(3, 4))

For non-contiguous replicate numbers, pass explicit replicate IDs:

.. code-block:: python

    study.add_process("drug", "CTRL", "DRUG", n_rep=((1, 3, 4), (1, 2, 5)))

Install FLiPPR
--------------

Use ``uv`` for a clean, reproducible Python environment:

.. code-block:: bash

    uv venv
    uv pip install flippr

For a source checkout of this repository:

.. code-block:: bash

    uv sync

DDA workflow
------------

Use this path for DDA LFQ LiP-MS data quantified by FragPipe and IonQuant.

1. Configure FragPipe v24.0
~~~~~~~~~~~~~~~~~~~~~~~~~~~

Download and launch FragPipe v24.0 from the `FragPipe releases`_ page. On
Windows, the installer bundles a Java runtime and Python. On Linux, FragPipe
v24.0 requires Java 11 or newer.

On the ``Config`` tab, confirm that FragPipe can find the tools it needs:
MSFragger, IonQuant, and the bundled or configured auxiliary tools. The official
FragPipe usage guide recommends starting from a built-in workflow and only then
customizing settings.

2. Choose a DDA LFQ workflow
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

On the ``Workflow`` tab, load a DDA label-free workflow appropriate for the
experiment, commonly ``LFQ-MBR`` or a closely related LFQ workflow.

For LiP-MS, the search settings matter. If the biological question depends on
semi-tryptic or non-tryptic cleavage products, make sure the MSFragger enzyme
and digestion settings can identify those peptides. Keep the LiP and TrP
workflows as similar as practical so normalization is interpretable.

Load all files that belong to the same FragPipe experiment together and give
the samples clear names:

.. code-block:: text

    CTRL_1.raw
    CTRL_2.raw
    CTRL_3.raw
    DRUG_1.raw
    DRUG_2.raw
    DRUG_3.raw

For Thermo files, the FragPipe documentation recommends converting ``.raw``
files to centroided ``.mzML`` with peak picking for many workflows. Direct raw
reading can work, but conversion makes runs easier to reproduce across systems.

3. Configure database and search settings
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

On the ``Database`` tab, select the FASTA database used for the experiment.
Use the same database for LiP and TrP runs. If you add contaminants or decoys,
keep that setup consistent across all output directories that will be compared
in FLiPPR.

On the MSFragger-related tabs, review:

- Enzyme and missed cleavage settings.
- Fixed and variable modifications.
- Precursor and fragment mass tolerances.
- Validation and FDR settings.

FLiPPR assumes the FragPipe reports have already been filtered and quantified.
It does not redo peptide-spectrum matching or protein inference.

4. Enable MS1 label-free quantification
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

On the ``Quant (MS1)`` tab, use IonQuant for label-free quantification.
FragPipe documents IonQuant as the default LFQ tool. For protein-level TrP
normalization in DDA, FLiPPR uses ``MaxLFQ Intensity`` by default from
``combined_protein.tsv``.

Recommended checks:

- Use MaxLFQ for protein quantification when available.
- Decide whether match-between-runs is appropriate for the design.
- Keep the MBR settings consistent between LiP and TrP analyses.
- Use the same sample naming pattern for every run.

5. Run FragPipe and inspect output
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Run FragPipe from the ``Run`` tab. After completion, the LiP output directory
should look like this:

.. code-block:: text

    LiP_DDA_FragPipe/
    ├── combined_ion.tsv
    ├── combined_protein.tsv
    ├── experiment_annotation.tsv
    ├── log_*.txt
    └── ...

If you have matched TrP data, run it into a separate output directory:

.. code-block:: text

    TrP_DDA_FragPipe/
    ├── combined_ion.tsv
    ├── combined_protein.tsv
    ├── experiment_annotation.tsv
    ├── log_*.txt
    └── ...

Check ``experiment_annotation.tsv`` before running FLiPPR. The sample names
should match the condition and replicate naming pattern that you will pass to
``add_process()``.

6. Run FLiPPR on DDA LiP data
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

    from pathlib import Path

    import flippr

    lip_dir = Path("LiP_DDA_FragPipe")

    study = flippr.Study(lip=lip_dir, method="dda")
    print(study.samples)

    study.add_process(
        pid="drug",
        lip_ctrl="CTRL",
        lip_test="DRUG",
        n_rep=3,
    )

    results = study.run()
    drug = results["drug"]

    drug.ion
    drug.modified_peptide
    drug.peptide
    drug.cut_site
    drug.protein_summary

7. Run FLiPPR with DDA TrP normalization
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

    from pathlib import Path

    import flippr

    lip_dir = Path("LiP_DDA_FragPipe")
    trp_dir = Path("TrP_DDA_FragPipe")

    study = flippr.Study(lip=lip_dir, trp=trp_dir, method="dda")
    print(study.samples)

    study.add_process(
        pid="drug_norm",
        lip_ctrl="CTRL",
        lip_test="DRUG",
        n_rep=3,
        trp_ctrl="CTRL",
        trp_test="DRUG",
        trp_n_rep=3,
    )

    results = study.run()
    drug = results["drug_norm"]

    drug.ion
    drug.trp_protein
    drug.protein_summary

DIA workflow
------------

Use this path for DIA data processed in FragPipe v24.0 with DIA-NN
quantification.

1. Choose a DIA workflow in FragPipe
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

FragPipe documents two main DIA workflows:

- ``DIA_SpecLib_Quant`` uses MSFragger-DIA to build a spectral library, then
  quantifies with DIA-NN.
- ``DIA_DIA-Umpire_SpecLib_Quant`` uses DIA-Umpire to generate pseudo-MS/MS
  spectra, searches those spectra in DDA mode, then quantifies with DIA-NN.

For most direct DIA searches, start with ``DIA_SpecLib_Quant``. If the data
were acquired with overlapping or staggered windows, follow the FragPipe DIA
documentation and demultiplex to ``.mzML`` before analysis.

If you already have a spectral library and only need quantification, use the
``Quant (DIA)`` tab to run only the DIA-NN quantification step with that
library.

2. Load files and verify data types
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

On the ``Workflow`` tab, load the DIA files and check the inferred data type.
FragPipe guesses from file and folder names, but you should verify it.

Use the same sample naming pattern used for DDA:

.. code-block:: text

    CTRL_1
    CTRL_2
    CTRL_3
    DRUG_1
    DRUG_2
    DRUG_3

If you include separate DDA files for spectral library construction, make sure
those files are marked with the correct DDA data type. FragPipe notes that
pseudo-MS/MS files from DIA-Umpire should be designated as DDA data.

3. Configure spectral library and DIA-NN quantification
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

If spectral library generation is enabled, confirm that Python and EasyPQP are
available in FragPipe. The FragPipe DIA tutorial notes that Python with EasyPQP
is required for spectral library generation.

On ``Quant (DIA)``, keep DIA-NN enabled. The generated or supplied library will
be passed to DIA-NN for quantification. Review DIA-NN settings such as library
source, quantification strategy, MBR-like options, and protein inference
options according to the experimental design.

4. Run FragPipe and inspect output
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

After a DIA run, FLiPPR expects this structure:

.. code-block:: text

    LiP_DIA_FragPipe/
    ├── ion.tsv
    ├── experiment_annotation.tsv
    ├── dia-quant-output/
    │   ├── report.pr_matrix.tsv
    │   └── report.pg_matrix.tsv
    ├── log_*.txt
    └── ...

The ``report.pr_matrix.tsv`` file contains precursor-level DIA-NN quantities.
FLiPPR combines those quantities with peptide/protein location metadata from
``ion.tsv``. The ``report.pg_matrix.tsv`` file is used for optional protein
normalization when a matched TrP DIA run is provided.

Some FragPipe/DIA-NN documentation and older workflows refer to a
``diann-output`` directory. FLiPPR currently looks for ``dia-quant-output``.
If your FragPipe v24 output uses a different folder name, keep the contents
unchanged and make the path match before loading it with FLiPPR.

5. Run FLiPPR on DIA LiP data
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

    from pathlib import Path

    import flippr

    lip_dir = Path("LiP_DIA_FragPipe")

    study = flippr.Study(lip=lip_dir, method="dia")
    print(study.samples)

    study.add_process(
        pid="drug",
        lip_ctrl="CTRL",
        lip_test="DRUG",
        n_rep=3,
    )

    results = study.run()
    drug = results["drug"]

    drug.ion
    drug.modified_peptide
    drug.peptide
    drug.cut_site
    drug.protein_summary

6. Run FLiPPR with DIA TrP normalization
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

    from pathlib import Path

    import flippr

    lip_dir = Path("LiP_DIA_FragPipe")
    trp_dir = Path("TrP_DIA_FragPipe")

    study = flippr.Study(lip=lip_dir, trp=trp_dir, method="dia")

    study.add_process(
        pid="drug_norm",
        lip_ctrl="CTRL",
        lip_test="DRUG",
        n_rep=3,
        trp_ctrl="CTRL",
        trp_test="DRUG",
        trp_n_rep=3,
    )

    results = study.run()
    drug = results["drug_norm"]

    drug.ion
    drug.trp_protein
    drug.protein_summary

Working with multiple comparisons
---------------------------------

A single ``Study`` can hold multiple comparisons as long as they use the same
FragPipe output directory.

.. code-block:: python

    study = flippr.Study(lip="LiP_DDA_FragPipe", trp="TrP_DDA_FragPipe", method="dda")

    study.add_process("low", "CTRL", "DRUG_LOW", n_rep=3, trp_ctrl="CTRL", trp_test="DRUG_LOW", trp_n_rep=3)
    study.add_process("high", "CTRL", "DRUG_HIGH", n_rep=3, trp_ctrl="CTRL", trp_test="DRUG_HIGH", trp_n_rep=3)

    results = study.run()

Exporting results
-----------------

FLiPPR returns Polars dataframes. Use any standard Polars writer:

.. code-block:: python

    result = results["drug"]

    result.ion.write_csv("drug_ion.csv")
    result.peptide.write_csv("drug_peptide.csv")
    result.cut_site.write_csv("drug_cut_site.csv")
    result.protein_summary.write_csv("drug_protein_summary.csv")

For larger studies, Parquet is usually faster and preserves types better:

.. code-block:: python

    result.ion.write_parquet("drug_ion.parquet")
    result.protein_summary.write_parquet("drug_protein_summary.parquet")

Important FLiPPR parameters
---------------------------

Global parameters live in ``flippr.rcParams``. Change them before calling
``study.run()``.

.. code-block:: python

    import flippr

    flippr.rcParams["ion.missing_intensity_thresh"] = 1
    flippr.rcParams["protein.fc_sig_thresh"] = 1.0
    flippr.rcParams["protein.pval_sig_thresh"] = 0.01
    flippr.rcParams["protein.adj_pval_sig_thresh"] = 0.05

Useful settings:

- ``ion.missing_intensity_thresh`` controls how many missing control/test
  replicate intensities are tolerated before an ion is filtered.
- ``ion.aon_impute_loc`` and ``ion.aon_impute_scale`` control all-or-none
  imputation values.
- ``trp_protein.intensity_value`` controls the DDA protein intensity column
  used for TrP normalization. The default is ``MaxLFQ Intensity``.
- ``protein.fc_sig_thresh``, ``protein.pval_sig_thresh``, and
  ``protein.adj_pval_sig_thresh`` control protein-summary significance counts.

Quality checks
--------------

Before trusting a run, check these items.

FragPipe output files:

- DDA output has ``combined_ion.tsv`` and ``combined_protein.tsv``.
- DIA output has ``ion.tsv`` and both DIA-NN matrices.
- ``experiment_annotation.tsv`` exists in every LiP and TrP output directory.
- Sample names match the ``CTRL_1`` / ``DRUG_1`` style expected by FLiPPR.

Replicates:

- ``n_rep=3`` means FLiPPR will look for ``CONDITION_1`` through
  ``CONDITION_3``.
- ``n_rep=(3, 4)`` means three control replicates and four test replicates.
- ``n_rep=((1, 3, 4), (1, 2, 5))`` means use explicit replicate labels.

Biology and design:

- LiP and TrP runs should use compatible FASTA, search, and quantification
  settings.
- TrP normalization should use matched biological conditions, not unrelated
  controls.
- MBR should be chosen with the experimental design in mind. For highly
  separated conditions, unrestricted transfer can add noise.

Troubleshooting
---------------

``FileNotFoundError`` when creating ``Study``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The output directory is missing one of the files FLiPPR expects. Check the file
tree against the DDA or DIA input lists above. For DIA, check the exact name of
the DIA-NN output directory.

``ColumnNotFoundError`` during ``study.run()``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The sample names or replicate numbers passed to ``add_process()`` do not match
the FragPipe intensity columns. Open ``experiment_annotation.tsv`` and the
header row of the relevant quantification table, then update ``lip_ctrl``,
``lip_test``, or ``n_rep``.

No rows after filtering
~~~~~~~~~~~~~~~~~~~~~~~

The missing-intensity filter may be too strict for the dataset, or the sample
names may point to the wrong columns. Check ``ion.missing_intensity_thresh`` and
verify the sample names.

Unexpected TrP normalization values
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For DDA, verify that ``combined_protein.tsv`` contains the intensity column
configured by ``flippr.rcParams["trp_protein.intensity_value"]``. For DIA,
FLiPPR uses DIA-NN protein-group intensities from ``report.pg_matrix.tsv``.

References
----------

- `FragPipe releases`_
- `FragPipe usage guide`_
- `FragPipe DIA tutorial`_
- `FragPipe output guide`_

.. _FragPipe releases: https://github.com/Nesvilab/FragPipe/releases
.. _FragPipe usage guide: https://fragpipe.nesvilab.org/docs/tutorial_fragpipe.html
.. _FragPipe DIA tutorial: https://fragpipe.nesvilab.org/docs/tutorial_DIA.html
.. _FragPipe output guide: https://fragpipe.nesvilab.org/docs/tutorial_fragpipe_outputs.html
