/* ==========================================
   PCA PLOT
========================================== */

export function drawPCA(data) {

    const container =
        document.getElementById(
            "pcaPlot"
        );

    if (!data || data.length === 0) {

        container.innerHTML =
            "<p>No PCA data available.</p>";

        return;
    }

    const groups =
        [...new Set(
            data.map(
                x => x.group
            )
        )];

    const traces = groups.map(group => {

        const subset =
            data.filter(
                d => d.group === group
            );

        return {

            x:
                subset.map(
                    d => d.PC1
                ),

            y:
                subset.map(
                    d => d.PC2
                ),

            text:
                subset.map(
                    d => d.sample
                ),

            customdata:
                subset,

            mode:
                "markers+text",

            textposition:
                "top center",

            type:
                "scatter",

            name:
                group,

            marker: {
                size: 14
            }
        };

    });

    const layout = {

        title:
            "Principal Component Analysis",

        xaxis: {
            title: "PC1"
        },

        yaxis: {
            title: "PC2"
        },

        hovermode:
            "closest"
    };

    Plotly.newPlot(
        "pcaPlot",
        traces,
        layout,
        {
            responsive: true
        }
    );
}

/* ==========================================
   VOLCANO PLOT
========================================== */

export function drawVolcano(data) {

    const container =
        document.getElementById(
            "volcanoPlot"
        );

    if (!data || data.length === 0) {

        container.innerHTML =
            "<p>No differential expression results.</p>";

        return;
    }

    const x =
        data.map(
            p => p.log2fc
        );

    const y =
        data.map(
            p =>
                -Math.log10(
                    Math.max(
                        p.adjP,
                        1e-300
                    )
                )
        );

    const colors =
        data.map(row => {

            if (
                row.significant &&
                row.log2fc > 1
            ) {
                return "#d62728";
            }

            if (
                row.significant &&
                row.log2fc < -1
            ) {
                return "#1f77b4";
            }

            return "#bdbdbd";
        });

    const trace = {

        x,

        y,

        mode:
            "markers",

        type:
            "scattergl",

        customdata:
            data,

        marker: {

            size: 8,

            color:
                colors,

            opacity:
                0.75
        },

        hovertemplate:
            "<b>%{customdata.gene}</b><br>" +
            "Accession: %{customdata.accession}<br>" +
            "Log2FC: %{x:.2f}<br>" +
            "Adj.P: %{customdata.adjP:.4e}<br>" +
            "<extra></extra>"
    };

    const layout = {

        title:
            "Volcano Plot",

        xaxis: {

            title:
                "Log2 Fold Change",

            zeroline:
                true
        },

        yaxis: {

            title:
                "-log10(FDR)"
        },

        shapes: [

            {
                type: "line",
                x0: -1,
                x1: -1,
                y0: 0,
                y1: Math.max(...y),
                line: {
                    dash: "dash",
                    color: "gray"
                }
            },

            {
                type: "line",
                x0: 1,
                x1: 1,
                y0: 0,
                y1: Math.max(...y),
                line: {
                    dash: "dash",
                    color: "gray"
                }
            },

            {
                type: "line",
                x0:
                    Math.min(...x),

                x1:
                    Math.max(...x),

                y0:
                    -Math.log10(0.05),

                y1:
                    -Math.log10(0.05),

                line: {
                    dash: "dash",
                    color: "gray"
                }
            }
        ],

        hovermode:
            "closest"
    };

    Plotly.newPlot(
        "volcanoPlot",
        [trace],
        layout,
        {
            responsive: true
        }
    );

    attachVolcanoClick();
}

/* ==========================================
   VOLCANO CLICK HANDLER
========================================== */

function attachVolcanoClick() {

    const volcano =
        document.getElementById(
            "volcanoPlot"
        );

    volcano.removeAllListeners?.(
        "plotly_click"
    );

    volcano.on(
        "plotly_click",

        event => {

            const protein =
                event.points[0]
                .customdata;

            showProteinInfo(
                protein
            );

        }
    );
}

/* ==========================================
   PROTEIN INFO PANEL
========================================== */

function showProteinInfo(
    protein
) {

    const panel =
        document.getElementById(
            "proteinInfo"
        );

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
            <th>Significant</th>
            <td>${protein.significant}</td>
        </tr>

    </table>

    `;
}

/* ==========================================
   HEATMAP
========================================== */

export function drawHeatmap(
    heatmapData
) {

    const container =
        document.getElementById(
            "heatmapPlot"
        );

    if (
        !heatmapData ||
        !heatmapData.zscores
    ) {

        container.innerHTML =
            "<p>No heatmap data available.</p>";

        return;
    }

    const trace = {

        z:
            heatmapData.zscores,

        x:
            heatmapData.samples,

        y:
            heatmapData.proteins,

        type:
            "heatmap",

        colorscale:
            "RdBu",

        reversescale:
            true,

        zmid:
            0,

        hovertemplate:
            "Protein: %{y}<br>" +
            "Sample: %{x}<br>" +
            "Z-score: %{z:.2f}<extra></extra>"
    };

    const layout = {

        title:
            "Protein Heatmap",

        height:
            900,

        yaxis: {

            automargin:
                true
        }
    };

    Plotly.newPlot(
        "heatmapPlot",
        [trace],
        layout,
        {
            responsive: true
        }
    );
}

/* ==========================================
   RESULTS TABLE
========================================== */

export function drawResultsTable(
    results
) {

    const table =
        document.getElementById(
            "resultsTable"
        );

    if (
        !results ||
        results.length === 0
    ) {

        table.innerHTML =
            "<tr><td>No results available.</td></tr>";

        return;
    }

    const header = `

    <tr>

        <th>Gene</th>
        <th>Accession</th>
        <th>Log2FC</th>
        <th>Adj.P</th>
        <th>Significant</th>

    </tr>

    `;

    const rows =
        results
        .slice(0,1000)
        .map(r => `

        <tr>

            <td>${r.gene}</td>

            <td>${r.accession}</td>

            <td>${r.log2fc.toFixed(3)}</td>

            <td>${r.adjP.toExponential(3)}</td>

            <td>${r.significant}</td>

        </tr>

        `)
        .join("");

    table.innerHTML =
        header + rows;
}

/* ==========================================
   DOWNLOAD PLOT PNG
========================================== */

export function downloadPlot(
    divId,
    filename
) {

    Plotly.downloadImage(
        divId,
        {
            format: "png",
            filename,
            width: 1600,
            height: 1200
        }
    );
}
