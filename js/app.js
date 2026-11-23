import { loadExcel } from "./analysis.js";

import {
    drawVolcano,
    drawPCA,
    drawHeatmap
} from "./plots.js";

import {
    buildBatchEditor,
    getBatchAssignments
} from "./ui.js";

import {
    runCombat
} from "./combat.js";

/* ---------------------------------------
   GLOBAL APPLICATION STATE
--------------------------------------- */

let analysisResults = null;

/* ---------------------------------------
   FILE UPLOAD
--------------------------------------- */

document
    .getElementById("fileInput")
    .addEventListener("change", async (e) => {

        const file = e.target.files[0];

        if (!file) return;

        try {

            analysisResults =
                await loadExcel(file);

            buildBatchEditor(
                analysisResults.metadata
            );

            document.getElementById(
                "combatStatus"
            ).innerHTML = "";

            render();

        } catch (error) {

            console.error(error);

            alert(
                "Failed to load Proteome Discoverer file."
            );
        }
    });

/* ---------------------------------------
   COMBAT BUTTON
--------------------------------------- */

document
    .getElementById("applyCombatBtn")
    .addEventListener("click", async () => {

        if (!analysisResults) {

            alert(
                "Please upload a PD file first."
            );

            return;
        }

        const btn =
            document.getElementById(
                "applyCombatBtn"
            );

        const status =
            document.getElementById(
                "combatStatus"
            );

        btn.disabled = true;

        status.innerHTML =
            "Running ComBat batch correction...";

        try {

            const batchAssignments =
                getBatchAssignments();

            const correctedMatrix =
                await runCombat(
                    analysisResults.logMatrix,
                    analysisResults.metadata,
                    batchAssignments
                );

            analysisResults.corrected = true;

            analysisResults.correctedMatrix =
                correctedMatrix;

            /*
             * Future:
             * Recompute downstream analyses
             * from corrected matrix
             */

            // analysisResults.correctedPCA =
            //     await runPCA(correctedMatrix);

            // analysisResults.correctedVolcano =
            //     await runDifferentialAnalysis(
            //         correctedMatrix,
            //         analysisResults.metadata
            //     );

            // analysisResults.correctedHeatmap =
            //     await runHeatmapAnalysis(
            //         correctedMatrix
            //     );

            status.innerHTML =
                "✓ ComBat completed successfully";

            render();

        } catch (error) {

            console.error(error);

            status.innerHTML =
                "✗ ComBat failed";

        } finally {

            btn.disabled = false;
        }
    });

/* ---------------------------------------
   RENDER
--------------------------------------- */

function render() {

    if (!analysisResults) return;

    /*
     * If corrected analyses exist,
     * display those instead of raw.
     */

    const pcaData =
        analysisResults.correctedPCA ||
        analysisResults.pca;

    const volcanoData =
        analysisResults.correctedVolcano ||
        analysisResults.volcano;

    const heatmapData =
        analysisResults.correctedHeatmap ||
        analysisResults.heatmap;

    if (pcaData) {
        drawPCA(pcaData);
    }

    if (volcanoData) {
        drawVolcano(volcanoData);
    }

    if (heatmapData) {
        drawHeatmap(heatmapData);
    }

    updateDataStatus();
}

/* ---------------------------------------
   DATA STATUS DISPLAY
--------------------------------------- */

function updateDataStatus() {

    const statusElement =
        document.getElementById(
            "dataStatus"
        );

    if (!statusElement) return;

    if (
        analysisResults &&
        analysisResults.corrected
    ) {

        statusElement.innerHTML =
            "Showing ComBat-corrected data";

    } else {

        statusElement.innerHTML =
            "Showing raw data";
    }
}
