import { loadExcel } from "./analysis.js";

import {
    drawVolcano,
    drawPCA,
    drawHeatmap,
    drawResultsTable,
    getCurrentVolcanoData
} from "./plots.js";

import {
    buildBatchEditor,
    getBatchAssignments,
    setupProteinSearch,
    displayProteinInfo,
    highlightProtein,
    updateComparisonLabel
} from "./ui.js";

import {
    runCombat
} from "./combat.js";

/* ==========================================
   GLOBAL STATE
========================================== */

let analysisResults = null;

/* ==========================================
   COMPARISON DROPDOWN
========================================== */

function populateComparisonDropdown() {

    const selector =
        document.getElementById(
            "comparisonSelector"
        );

    if (!selector) return;

    selector.innerHTML = "";

    const comparisons =
        Object.keys(
            analysisResults.volcano
        );

    comparisons.forEach(name => {

        const option =
            document.createElement(
                "option"
            );

        option.value = name;
        option.textContent = name;

        selector.appendChild(
            option
        );

    });

}

/* ==========================================
   COMPARISON CHANGE
========================================== */

document
    .getElementById(
        "comparisonSelector"
    )
    .addEventListener(
        "change",
        () => {

            const comparison =
                document.getElementById(
                    "comparisonSelector"
                ).value;

            updateComparisonLabel(
                comparison
            );

            render();

        }
    );

/* ==========================================
   FILE UPLOAD
========================================== */

document
    .getElementById("fileInput")
    .addEventListener(
        "change",
        async (e) => {

            const file =
                e.target.files[0];

            if (!file) return;

            try {

                analysisResults =
                    await loadExcel(
                        file
                    );

                buildBatchEditor(
                    analysisResults.metadata
                );

                populateComparisonDropdown();

                const selector =
                    document.getElementById(
                        "comparisonSelector"
                    );

                const firstComparison =
                    Object.keys(
                        analysisResults.volcano
                    )[0];

                if (
                    selector &&
                    firstComparison
                ) {

                    selector.value =
                        firstComparison;

                    updateComparisonLabel(
                        firstComparison
                    );
                }

                setupProteinSearch(

                    getCurrentVolcanoData,

                    protein => {

                        displayProteinInfo(
                            protein
                        );

                        highlightProtein(
                            protein.accession
                        );

                    }

                );

                document.getElementById(
                    "combatStatus"
                ).textContent =
                    "Raw data loaded";

                render();

            } catch(error) {

                console.error(
                    error
                );

                alert(
                    "Failed to load Proteome Discoverer file."
                );

            }

        }
    );

/* ==========================================
   COMBAT
========================================== */

document
    .getElementById(
        "applyCombatBtn"
    )
    .addEventListener(
        "click",

        async () => {

            if (
                !analysisResults
            ) {

                alert(
                    "Please load a file first."
                );

                return;

            }

            const button =
                document.getElementById(
                    "applyCombatBtn"
                );

            const status =
                document.getElementById(
                    "combatStatus"
                );

            button.disabled = true;

            status.textContent =
                "Running ComBat batch correction...";

            try {

                const batches =
                    getBatchAssignments();

                const correctedMatrix =
                    await runCombat(

                        analysisResults.logMatrix,

                        analysisResults.metadata,

                        batches

                    );

                analysisResults.corrected =
                    true;

                analysisResults.correctedMatrix =
                    correctedMatrix;

                /*
                 * Future enhancement:
                 *
                 * Recalculate PCA
                 * Recalculate DE
                 * Recalculate Volcano
                 * Recalculate Heatmap
                 *
                 * from corrected matrix
                 */

                status.textContent =
                    "✓ ComBat completed successfully";

                render();

            } catch(error) {

                console.error(
                    error
                );

                status.textContent =
                    "✗ ComBat failed";

            } finally {

                button.disabled =
                    false;

            }

        }
    );

/* ==========================================
   RENDER
========================================== */

function render() {

    if (
        !analysisResults
    ) {
        return;
    }

    const selector =
        document.getElementById(
            "comparisonSelector"
        );

    const comparison =
        selector?.value ||
        Object.keys(
            analysisResults.volcano
        )[0];

    if (
        comparison
    ) {

        updateComparisonLabel(
            comparison
        );

    }

    const volcanoData =
        analysisResults.volcano[
            comparison
        ];

    const resultsData =
        analysisResults.results[
            comparison
        ];

    const pcaData =
        analysisResults.correctedPCA ||
        analysisResults.pca;

    const heatmapData =
        analysisResults.correctedHeatmap ||
        analysisResults.heatmap;

    if (
        pcaData
    ) {

        drawPCA(
            pcaData
        );

    }

    if (
        volcanoData
    ) {

        drawVolcano(
            volcanoData
        );

    }

    if (
        heatmapData
    ) {

        drawHeatmap(
            heatmapData
        );

    }

    if (
        resultsData
    ) {

        drawResultsTable(
            resultsData
        );

    }

    updateDataStatus();

}

/* ==========================================
   DATA STATUS
========================================== */

function updateDataStatus() {

    const status =
        document.getElementById(
            "dataStatus"
        );

    if (!status) return;

    if (
        analysisResults &&
        analysisResults.corrected
    ) {

        status.textContent =
            "Showing ComBat-corrected data";

    } else {

        status.textContent =
            "Showing raw data";

    }

}
