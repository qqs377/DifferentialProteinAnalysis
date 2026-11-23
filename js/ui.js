/* ==========================================
   BATCH EDITOR
========================================== */

export function buildBatchEditor(metadata) {

    const container =
        document.getElementById(
            "batchEditor"
        );

    if (!container) return;

    container.innerHTML = "";

    const table =
        document.createElement(
            "table"
        );

    table.className =
        "batch-table";

    table.innerHTML = `

        <thead>

            <tr>

                <th>Sample</th>
                <th>Group</th>
                <th>Batch</th>

            </tr>

        </thead>

        <tbody>

        ${metadata.map(item => `

            <tr>

                <td>${item.sample}</td>

                <td>${item.group}</td>

                <td>

                    <select
                        class="batch-select"
                        data-sample="${item.sample}"
                    >

                        <option value="Batch1">
                            Batch1
                        </option>

                        <option value="Batch2">
                            Batch2
                        </option>

                        <option value="Batch3">
                            Batch3
                        </option>

                        <option value="Batch4">
                            Batch4
                        </option>

                    </select>

                </td>

            </tr>

        `).join("")}

        </tbody>

    `;

    container.appendChild(table);

    addBatchCSVUploader();

    const btn =
        document.getElementById(
            "applyCombatBtn"
        );

    if (btn) {
        btn.disabled = false;
    }
}

/* ==========================================
   GET BATCH ASSIGNMENTS
========================================== */

export function getBatchAssignments() {

    const assignments = {};

    document
        .querySelectorAll(
            ".batch-select"
        )
        .forEach(select => {

            assignments[
                select.dataset.sample
            ] = select.value;

        });

    return assignments;
}

/* ==========================================
   CSV IMPORTER
========================================== */

function addBatchCSVUploader() {

    const container =
        document.getElementById(
            "batchEditor"
        );

    if (
        document.getElementById(
            "batchCsvInput"
        )
    ) {
        return;
    }

    const uploadDiv =
        document.createElement(
            "div"
        );

    uploadDiv.innerHTML = `

        <br>

        <label>

            Import Batch Metadata

        </label>

        <br>

        <input
            type="file"
            id="batchCsvInput"
            accept=".csv"
        >

    `;

    container.appendChild(
        uploadDiv
    );

    document
        .getElementById(
            "batchCsvInput"
        )
        .addEventListener(
            "change",
            loadBatchCSV
        );
}

/* ==========================================
   LOAD BATCH CSV
========================================== */

function loadBatchCSV(event) {

    const file =
        event.target.files[0];

    if (!file) return;

    const reader =
        new FileReader();

    reader.onload =
        function(e) {

            const text =
                e.target.result;

            const rows =
                text.trim()
                    .split("\n");

            rows
                .slice(1)
                .forEach(row => {

                    const parts =
                        row.split(",");

                    if (
                        parts.length < 2
                    ) {
                        return;
                    }

                    const sample =
                        parts[0]
                        .trim();

                    const batch =
                        parts[1]
                        .trim();

                    const selector =
                        document.querySelector(

                            `.batch-select[data-sample="${sample}"]`

                        );

                    if (
                        selector
                    ) {

                        selector.value =
                            batch;

                    }

                });

        };

    reader.readAsText(
        file
    );
}

/* ==========================================
   PROTEIN SEARCH
========================================== */

export function setupProteinSearch(
    getCurrentVolcanoData,
    onProteinFound
) {

    const button =
        document.getElementById(
            "searchProteinBtn"
        );

    if (!button) return;

    button.onclick =
        () => {

            const query =
                document
                .getElementById(
                    "proteinSearch"
                )
                .value
                .trim()
                .toUpperCase();

            if (
                !query
            ) {
                return;
            }

            const volcanoData =
                getCurrentVolcanoData();

            if (
                !volcanoData
            ) {
                return;
            }

            const match =
                volcanoData.find(
                    protein =>

                        protein.gene
                            ?.toUpperCase()
                            .includes(
                                query
                            )

                        ||

                        protein.accession
                            ?.toUpperCase()
                            .includes(
                                query
                            )
                );

            if (
                !match
            ) {

                alert(
                    "Protein not found."
                );

                return;

            }

            if (
                onProteinFound
            ) {

                onProteinFound(
                    match
                );

            }

        };

}

/* ==========================================
   HIGHLIGHT PROTEIN
========================================== */

export function highlightProtein(
    accession
) {

    const plot =
        document.getElementById(
            "volcanoPlot"
        );

    if (
        !plot ||
        !plot.data
    ) {
        return;
    }

    const trace =
        plot.data[0];

    if (
        !trace.customdata
    ) {
        return;
    }

    const sizes =
        trace.customdata.map(
            protein =>

            protein.accession ===
            accession

            ? 18
            : 8
        );

    const colors =
        trace.customdata.map(
            protein =>

            protein.accession ===
            accession

            ? "#FFD700"
            : trace.marker.color[
                trace.customdata.indexOf(
                    protein
                )
            ]
        );

    Plotly.restyle(
        plot,
        {
            "marker.size":[
                sizes
            ],
            "marker.color":[
                colors
            ]
        }
    );
}

/* ==========================================
   COMPARISON LABEL
========================================== */

export function updateComparisonLabel(
    comparison
) {

    let label =
        document.getElementById(
            "comparisonLabel"
        );

    if (!label) {

        label =
            document.createElement(
                "div"
            );

        label.id =
            "comparisonLabel";

        const selector =
            document.getElementById(
                "comparisonSelector"
            );

        if (
            selector &&
            selector.parentNode
        ) {

            selector.parentNode.insertBefore(
                label,
                selector.nextSibling
            );

        }

    }

    label.innerHTML =

        `<b>Current Comparison:</b> ${comparison}`;

}

/* ==========================================
   PROTEIN INFO PANEL
========================================== */

export function displayProteinInfo(
    protein
) {

    const panel =
        document.getElementById(
            "proteinInfo"
        );

    if (!panel) return;

    panel.innerHTML = `

        <h3>

            ${protein.gene || "Unknown"}

        </h3>

        <table>

            <tr>
                <th>Accession</th>
                <td>${protein.accession}</td>
            </tr>

            <tr>
                <th>Group 1</th>
                <td>${protein.group1}</td>
            </tr>

            <tr>
                <th>Group 2</th>
                <td>${protein.group2}</td>
            </tr>

            <tr>
                <th>Log2FC</th>
                <td>${protein.log2fc.toFixed(3)}</td>
            </tr>

            <tr>
                <th>P-value</th>
                <td>${protein.pvalue.toExponential(3)}</td>
            </tr>

            <tr>
                <th>Adjusted P</th>
                <td>${protein.adjP.toExponential(3)}</td>
            </tr>

            <tr>
                <th>Status</th>
                <td>${protein.status}</td>
            </tr>

        </table>

    `;
}

/* ==========================================
   DOWNLOAD CSV
========================================== */

export function downloadCSV(
    rows,
    filename
) {

    if (
        !rows ||
        rows.length === 0
    ) {
        return;
    }

    const headers =
        Object.keys(
            rows[0]
        );

    const csv = [

        headers.join(","),

        ...rows.map(
            row =>

            headers.map(
                h =>
                JSON.stringify(
                    row[h]
                )
            ).join(",")
        )

    ].join("\n");

    const blob =
        new Blob(
            [csv],
            {
                type:
                    "text/csv;charset=utf-8;"
            }
        );

    const url =
        URL.createObjectURL(
            blob
        );

    const link =
        document.createElement(
            "a"
        );

    link.href = url;

    link.download =
        filename;

    document.body.appendChild(
        link
    );

    link.click();

    document.body.removeChild(
        link
    );

    URL.revokeObjectURL(
        url
    );
}
