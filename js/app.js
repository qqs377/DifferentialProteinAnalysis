import {
loadExcel
}
from "./analysis.js";

import {
drawVolcano,
drawPCA,
drawHeatmap
}
from "./plots.js";

let analysisResults;

document
.getElementById("fileInput")
.addEventListener(
"change",
async(e)=>{

const file=e.target.files[0];

analysisResults =
await loadExcel(file);

render();

}
);

async function render(){

drawPCA(
analysisResults.pca
);

drawVolcano(
analysisResults.volcano
);

drawHeatmap(
analysisResults.heatmap
);

}
