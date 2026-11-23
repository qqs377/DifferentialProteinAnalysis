import {
initializePyodide
}
from "./pyodide.js";

export async function runCombat(
abundanceMatrix,
metadata,
batches
){

const py =
await initializePyodide();

py.globals.set(
"abundance_data",
abundanceMatrix
);

py.globals.set(
"groups",
metadata.map(
m=>m.group
)
);

py.globals.set(
"batch_labels",
metadata.map(
m=>batches[m.sample]
)
);

const result =
await py.runPythonAsync(`

import pandas as pd
import numpy as np

try:

    from neuroCombat import neuroCombat

except:

    raise Exception(
        "neuroCombat not installed"
    )

df = pd.DataFrame(
    abundance_data
)

combat_result =
neuroCombat(
    dat=df.values,
    covars=pd.DataFrame({

        "Batch": batch_labels,
        "Group": groups

    }),
    batch_col="Batch",
    categorical_cols=["Group"]

)

combat_result["data"].tolist()

`);

return result;

}
