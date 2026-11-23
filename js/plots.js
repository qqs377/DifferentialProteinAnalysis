import {
    displayProteinInfo,
    highlightProtein
} from "./ui.js";

/* ==========================================
   CURRENT VOLCANO DATA
========================================== */

let currentVolcanoData = [];

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

    const groups = [
        ...new Set(
            data.map(
                d => d.group
            )
        )
    ];

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

    Plotly.newPlot(

        "pcaPlot",

        traces,

        {

            title:
                "Principal Component Analysis",

            xaxis: {

                title:
                    "PC1"

            },

            yaxis: {

                title:
                    "PC2"

            },

            hovermode:
                "closest"

        },

        {

            responsive:
                true

        }

    );

}

/* ==========================================
   VOLCANO
========================================== */

export function drawVolcano(data) {

    currentVolcanoData =
        data || [];

    const volcano =
        document.getElementById(
            "volcanoPlot"
        );

    if (
        !data ||
        data.length === 0
    ) {

        volcano.innerHTML =
            "<p>No volcano results available.</p>";

        return;

    }

    const colors =
        data.map(row => {

            if (
                row.status ===
                "Upregulated"
            )
                return "#d62728";

            if (
                row.status ===
                "Downregulated"
            )
                return "#1f77b4";

            return "#bdbdbd";

        });

    const trace = {

        x:
            data.map(
                x => x.log2fc
            ),

        y:
            data.map(
                x =>
                x.negLog10AdjP
            ),

        mode:
            "markers",

        type:
            "scattergl",

        customdata:
            data,

        marker: {

            color:
                colors,

            size: 8,

            opacity: 0.8

        },

        hovertemplate:

            "<b>%{customdata.gene}</b><br>" +

            "Accession: %{customdata.accession}<br>" +

            "Log2FC: %{customdata.log2fc:.3f}<br>" +

            "Adj.P: %{customdata.adjP:.2e}<br>" +

            "<extra></extra>"

    };

    const maxY =
        Math.max(
            ...data.map(
                x =>
                x.negLog10AdjP
            )
        );

    Plotly.newPlot(

        volcano,

        [trace],

        {

            title:
                "Volcano Plot",

            xaxis: {

                title:
                    "Log2 Fold Change"

            },

            yaxis: {

                title:
                    "-Log10 FDR"

            },

            shapes: [

                {
                    type: "line",

                    x0: -1,
                    x1: -1,

                    y0: 0,
                    y1: maxY,

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
                    y1: maxY,

                    line: {
                        dash: "dash",
                        color: "gray"
                    }
                },

                {
                    type: "line",

                    x0:
                        Math.min(
                            ...data.map(
                                x =>
                                x.log2fc
                            )
                        ),

                    x1:
                        Math.max(
                            ...data.map(
                                x =>
                                x.log2fc
                            )
                        ),

                    y0:
                        -Math.log10(
                            0.05
                        ),

                    y1:
                        -Math.log10(
                            0.05
                        ),

                    line: {

                        dash:
                            "dash",

                        color:
                            "gray"

                    }

                }

            ]

        },

        {

            responsive:
                true

        }

    );

    volcano.on(
        "plotly_click",
        event => {

            const protein =
                event.points[0]
                    .customdata;

            displayProteinInfo(
                protein
            );

            highlightProtein(
                protein.accession
            );

        }
    );

}

/* ==========================================
   PROVIDE DATA FOR SEARCH
========================================== */

export function getCurrentVolcanoData() {

    return currentVolcanoData;

}

/* ==========================================
   HEATMAP
========================================== */

export function drawHeatmap(
    heatmapData
) {

    const heatmap =
        document.getElementById(
            "heatmapPlot"
        );

    if (
        !heatmapData ||
        !heatmapData.zscores
    ) {

        heatmap.innerHTML =
            "<p>No heatmap data available.</p>";

        return;

    }

    Plotly.newPlot(

        heatmap,

        [

            {

                z:
                    heatmapData
                        .zscores,

                x:
                    heatmapData
                        .samples,

                y:
                    heatmapData
                        .proteins,

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

            }

        ],

        {

            title:
                "Protein Heatmap",

            height:
                900

        },

        {

            responsive:
                true

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

    table.innerHTML = `

    <thead>

        <tr>

            <th>Gene</th>

            <th>Accession</th>

            <th>Log2FC</th>

            <th>Adj.P</th>

            <th>Status</th>

        </tr>

    </thead>

    <tbody>

        ${results
            .slice(0,1000)
            .map(row => `

                <tr>

                    <td>${row.gene}</td>

                    <td>${row.accession}</td>

                    <td>${row.log2fc.toFixed(3)}</td>

                    <td>${row.adjP.toExponential(3)}</td>

                    <td>${row.status}</td>

                </tr>

            `)
            .join("")}

    </tbody>

    `;

}

/* ==========================================
   DOWNLOAD PLOT
========================================== */

export function downloadPlot(
    divId,
    filename
) {

    Plotly.downloadImage(

        divId,

        {

            format:
                "png",

            filename,

            width:
                1600,

            height:
                1200

        }

    );

}
