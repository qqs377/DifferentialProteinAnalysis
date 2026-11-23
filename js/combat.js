import {
    initializePyodide
} from "./pyodide.js";

import {
    runPCA,
    runDifferentialAnalysis,
    prepareHeatmap
} from "./analysis.js";

/* ==========================================
   MAIN COMBAT WORKFLOW
========================================== */

export async function runCombatWorkflow(
    analysisResults,
    batchAssignments
) {

    const correctedMatrix =
        await runCombat(

            analysisResults.logMatrix,

            analysisResults.metadata,

            batchAssignments

        );

    const correctedPCA =
        runPCA(
            correctedMatrix,
            analysisResults.metadata
        );

    const correctedDE =
        runDifferentialAnalysis(
            correctedMatrix,
            analysisResults.metadata
        );

    const correctedHeatmap =
        prepareHeatmap(
            correctedMatrix
        );

    return {

        correctedMatrix,

        correctedPCA,

        correctedVolcano:
            correctedDE.volcano,

        correctedResults:
            correctedDE.results,

        correctedHeatmap

    };

}

/* ==========================================
   COMBAT CORRECTION
========================================== */

export async function runCombat(
    matrix,
    metadata,
    batchAssignments
) {

    const pyodide =
        await initializePyodide();

    pyodide.globals.set(
        "batchAssignments",
        batchAssignments
    );

    pyodide.globals.set(
        "proteinMatrix",
        JSON.stringify(matrix)
    );

    pyodide.globals.set(
        "metadataRows",
        JSON.stringify(metadata)
    );

    const corrected =
        await pyodide.runPythonAsync(`

import json
import numpy as np
import pandas as pd

proteins = json.loads(
    proteinMatrix
)

metadata = json.loads(
    metadataRows
)

# -----------------------------------
# Build dataframe
# -----------------------------------

sample_columns = list(
    proteins[0]["values"].keys()
)

matrix_rows = []

for protein in proteins:

    row = []

    for sample in sample_columns:

        value = protein["values"].get(sample)

        if value is None:

            row.append(np.nan)

        else:

            row.append(float(value))

    matrix_rows.append(row)

df = pd.DataFrame(
    matrix_rows,
    columns=sample_columns
)

# -----------------------------------
# Metadata
# -----------------------------------

sample_to_batch = {}

for sample, batch in batchAssignments.items():

    sample_to_batch[sample] = batch

column_batches = []

for m in metadata:

    sample_name = m["sample"]

    column_batches.append(
        sample_to_batch.get(
            sample_name,
            "Batch1"
        )
    )

# -----------------------------------
# Median imputation
# -----------------------------------

df = df.T

df = df.fillna(
    df.median()
)

# -----------------------------------
# Simple batch centering
# -----------------------------------
#
# Replace with neuroCombat later.
#
# Subtract batch mean
# Add global mean
#
# -----------------------------------

global_mean = df.mean()

corrected = df.copy()

for batch in set(column_batches):

    sample_mask = [

        b == batch
        for b in column_batches

    ]

    batch_samples = corrected.index[
        sample_mask
    ]

    batch_mean = corrected.loc[
        batch_samples
    ].mean()

    corrected.loc[
        batch_samples
    ] = (

        corrected.loc[
            batch_samples
        ]

        - batch_mean

        + global_mean

    )

corrected = corrected.T

result = []

for i, protein in enumerate(proteins):

    values = {}

    for j, sample in enumerate(sample_columns):

        value = corrected.iloc[
            i,
            j
        ]

        if pd.isna(value):

            values[sample] = None

        else:

            values[sample] = float(value)

    result.append({

        "accession":
            protein["accession"],

        "gene":
            protein["gene"],

        "values":
            values

    })

json.dumps(result)

`);

    return JSON.parse(
        corrected
    );

}
