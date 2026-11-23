export function drawVolcano(df){

const trace={

x:df.map(p=>p.log2fc),

y:df.map(
p=>-Math.log10(p.adjP)
),

mode:"markers",

type:"scatter",

customdata:df,

marker:{

size:8

}

};

Plotly.newPlot(
"volcanoPlot",
[trace]
);

document
.getElementById("volcanoPlot")
.on(
"plotly_click",
function(event){

const protein =
event.points[0]
.customdata;

document
.getElementById(
"proteinInfo"
)
.innerHTML =

`
<h3>${protein.gene}</h3>

<p>
<b>Accession:</b>
${protein.accession}
</p>

<p>
<b>Log2FC:</b>
${protein.log2fc.toFixed(3)}
</p>

<p>
<b>Adjusted P:</b>
${protein.adjP}
</p>
`;

}
);

}
