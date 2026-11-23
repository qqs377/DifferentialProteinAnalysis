import PCA from "https://cdn.skypack.dev/ml-pca";

import jstat from "https://esm.sh/jstat";
const { jStat } = jstat;

/* ==========================================
   MAIN ENTRY POINT
========================================== */

export async function loadExcel(file) {

    const data =
        await file.arrayBuffer();

    const workbook =
        XLSX.read(data, {
            type: "array"
        });

    const sheet =
        workbook.Sheets[
            workbook.SheetNames[0]
        ];

    const rawRows =
        XLSX.utils.sheet_to_json(sheet);

    const abundanceCols =
        extractAbundanceColumns(rawRows);

    const metadata =
        extractSampleMetadata(
            abundanceCols
        );

    const filtered =
        filterProteins(rawRows);

    const abundanceMatrix =
        createAbundanceMatrix(
            filtered,
            abundanceCols
        );

    const logMatrix =
        logTransform(
            abundanceMatrix
        );

    const pca =
        runPCA(
            logMatrix,
            metadata
        );

    const differential =
        runDifferentialAnalysis(
            logMatrix,
            metadata
        );

    const heatmap =
        prepareHeatmap(
            logMatrix
        );

    return {

        rawRows,

        metadata,

        abundanceCols,

        logMatrix,

        pca,

        volcano:
            differential.volcano,

        results:
            differential.results,

        heatmap,

        corrected: false
    };
}

/* ==========================================
   ABUNDANCE COLUMNS
========================================== */

function extractAbundanceColumns(rows) {

    if (!rows.length)
        return [];

    return Object.keys(rows[0])
        .filter(col =>
            col.includes(
                "Abundances (Normalized)"
            )
        );
}

/* ==========================================
   SAMPLE METADATA
========================================== */

function extractSampleMetadata(
    abundanceCols
) {

    return abundanceCols.map(col => {

        const match =
            col.match(
                /:\s*([^:]+):\s*Sample,\s*(.+)$/
            );

        return {

            column: col,

            sample:
                match ?
                match[1].trim()
                : col,

            group:
                match ?
                match[2].trim()
                : "Unknown"
        };

    });
}

/* ==========================================
   FILTER PROTEINS
========================================== */

function filterProteins(rows) {

    return rows.filter(row => {

        const master =
            String(
                row["Master"] || ""
            ).toLowerCase();

        const contaminant =
            String(
                row["Contaminant"] || ""
            ).toLowerCase();

        const confidence =
            String(
                row["Protein FDR Confidence: Combined"] || ""
            ).toLowerCase();

        const peptides =
            Number(
                row["# Unique Peptides"]
            );

        if (
            master !==
            "master protein"
        ) {
            return false;
        }

        if (
            ["true", "1", "yes"]
            .includes(contaminant)
        ) {
            return false;
        }

        if (
            confidence !== "high"
        ) {
            return false;
        }

        if (
            isNaN(peptides)
        ) {
            return false;
        }

        if (
            peptides < 2
        ) {
            return false;
        }

        return true;
    });
}

/* ==========================================
   ABUNDANCE MATRIX
========================================== */

function createAbundanceMatrix(
    rows,
    abundanceCols
) {

    const matrix = [];

    rows.forEach(row => {

        const accession =
            row["Accession"];

        const gene =
            row["Gene Symbol"] ||
            "";

        const values = {};

        abundanceCols.forEach(col => {

            let v =
                Number(row[col]);

            if (
                isNaN(v) ||
                v <= 0
            ) {
                v = null;
            }

            values[col] = v;
        });

        matrix.push({

            accession,

            gene,

            values
        });

    });

    return matrix;
}

/* ==========================================
   LOG2 TRANSFORM
========================================== */

function logTransform(
    abundanceMatrix
) {

    return abundanceMatrix.map(
        protein => {

            const transformed =
                {};

            Object.keys(
                protein.values
            ).forEach(col => {

                const value =
                    protein.values[col];

                transformed[col] =
                    value
                        ? Math.log2(value)
                        : null;

            });

            return {

                accession:
                    protein.accession,

                gene:
                    protein.gene,

                values:
                    transformed
            };
        }
    );
}

/* ==========================================
   MEDIAN IMPUTATION
========================================== */

function median(arr) {

    const sorted =
        [...arr].sort(
            (a,b)=>a-b
        );

    const mid =
        Math.floor(
            sorted.length / 2
        );

    return sorted.length % 2
        ? sorted[mid]
        : (
            sorted[mid-1]
            +
            sorted[mid]
          ) / 2;
}

function imputeMissing(matrix) {

    const columns =
        Object.keys(
            matrix[0].values
        );

    const medians = {};

    columns.forEach(col => {

        const vals =
            matrix
            .map(
                p => p.values[col]
            )
            .filter(
                v => v !== null
            );

        medians[col] =
            median(vals);
    });

    return matrix.map(p => {

        const vals = {};

        columns.forEach(col => {

            vals[col] =
                p.values[col] ??
                medians[col];

        });

        return {

            accession:
                p.accession,

            gene:
                p.gene,

            values:
                vals
        };
    });
}

/* ==========================================
   PCA
========================================== */

export function runPCA(
    matrix,
    metadata
) {

    const imputed =
        imputeMissing(matrix);

    const columns =
        Object.keys(
            imputed[0].values
        );

    const sampleMap = {};

    columns.forEach(col => {

        sampleMap[col] = [];

        imputed.forEach(p => {

            sampleMap[col].push(
                p.values[col]
            );

        });

    });

    const X =
        Object.values(
            sampleMap
        );

    const pca =
        new PCA(X);

    const scores =
        pca.predict(X);

    const result = [];

    columns.forEach(
        (sample,i)=>{

        const meta =
            metadata.find(
                m =>
                m.column === sample
            );

        result.push({

            sample:
                meta.sample,

            group:
                meta.group,

            PC1:
                scores.data[i][0],

            PC2:
                scores.data[i][1]

        });

    });

    return result;
}

/* ==========================================
   DIFFERENTIAL EXPRESSION
========================================== */

export function runDifferentialAnalysis(
    matrix,
    metadata
) {

    const groups = [
        ...new Set(
            metadata.map(
                m => m.group
            )
        )
    ];

    if (groups.length < 2) {

        return {
            volcano: {},
            results: {}
        };

    }

    const comparisons = {};

    for (let i = 0; i < groups.length; i++) {

        for (let j = i + 1; j < groups.length; j++) {

            const g1 = groups[i];
            const g2 = groups[j];

            const cols1 =
                metadata
                .filter(
                    m => m.group === g1
                )
                .map(
                    m => m.column
                );

            const cols2 =
                metadata
                .filter(
                    m => m.group === g2
                )
                .map(
                    m => m.column
                );

            const results = [];

            matrix.forEach(protein => {

                const values1 =
                    cols1
                    .map(
                        c =>
                        protein.values[c]
                    )
                    .filter(
                        v =>
                        v !== null &&
                        !isNaN(v)
                    );

                const values2 =
                    cols2
                    .map(
                        c =>
                        protein.values[c]
                    )
                    .filter(
                        v =>
                        v !== null &&
                        !isNaN(v)
                    );

                if (
                    values1.length < 2 ||
                    values2.length < 2
                ) {
                    return;
                }

                const pvalue =
                    welchsTTest(
                        values1,
                        values2
                    );

                const mean1 =
                    jStat.mean(values1);

                const mean2 =
                    jStat.mean(values2);

                const log2fc =
                    mean1 - mean2;

                results.push({

                    accession:
                        protein.accession,

                    gene:
                        protein.gene || "",

                    group1:
                        g1,

                    group2:
                        g2,

                    nGroup1:
                        values1.length,

                    nGroup2:
                        values2.length,

                    meanGroup1:
                        mean1,

                    meanGroup2:
                        mean2,

                    log2fc,

                    foldChange:
                        Math.pow(
                            2,
                            log2fc
                        ),

                    pvalue
                });

            });

            applyBH(results);

            results.forEach(r => {

                r.significant =
                    r.adjP < 0.05 &&
                    Math.abs(
                        r.log2fc
                    ) >= 1;

                if (
                    r.significant &&
                    r.log2fc > 0
                ) {

                    r.status =
                        "Upregulated";

                } else if (
                    r.significant &&
                    r.log2fc < 0
                ) {

                    r.status =
                        "Downregulated";

                } else {

                    r.status =
                        "Not Significant";

                }

                r.negLog10AdjP =
                    -Math.log10(
                        Math.max(
                            r.adjP,
                            1e-300
                        )
                    );

            });

            results.sort(
                (a, b) =>
                a.adjP - b.adjP
            );

            const comparisonName =
                `${g1} vs ${g2}`;

            comparisons[
                comparisonName
            ] = results;

        }

    }

    return {

        volcano:
            comparisons,

        results:
            comparisons

    };

}

/* ==========================================
   WELCH TEST
========================================== */

function welchsTTest(a,b){

    const t =
        jStat.ttest(
            jStat.mean(a)
            -
            jStat.mean(b),

            a.length
            +
            b.length
            -
            2,

            2
        );

    return Math.abs(t);
}

/* ==========================================
   BH FDR
========================================== */

function applyBH(results){

    const pvals =
        results.map(
            r =>
            r.pvalue
        );

    const n =
        pvals.length;

    const sorted =
        [...pvals]
        .map(
          (p,i)=>({p,i})
        )
        .sort(
          (a,b)=>a.p-b.p
        );

    sorted.forEach(
        (item,rank)=>{

        const adj =
            (item.p*n)
            /
            (rank+1);

        results[
            item.i
        ].adjP =
            Math.min(
                adj,
                1
            );

    });
}

/* ==========================================
   HEATMAP
========================================== */

export function prepareHeatmap(
    matrix
) {

    const variances =
        matrix.map(p => {

        const vals =
            Object.values(
                p.values
            )
            .filter(
                v => v !== null
            );

        return {

            protein:p,

            variance:
                jStat.variance(
                    vals
                )
        };

    });

    variances.sort(
        (a,b)=>
        b.variance
        -
        a.variance
    );

    const top =
        variances
        .slice(0,50);

    const proteins =
        top.map(
            x =>
            x.protein.gene
            ||
            x.protein.accession
        );

    const samples =
        Object.keys(
            top[0]
            .protein
            .values
        );

    const zscores =
        top.map(item => {

            const vals =
                Object.values(
                    item.protein.values
                );

            const mean =
                jStat.mean(
                    vals
                );

            const sd =
                jStat.stdev(
                    vals
                );

            return vals.map(v => {

                if(sd === 0)
                    return 0;

                return (
                    v - mean
                ) / sd;

            });

        });

    return {

        proteins,

        samples,

        zscores

    };
}
