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

import {
buildBatchEditor
}
from "./ui.js";

import {
runCombat
}
from "./combat.js";

import {
getBatchAssignments
}
from "./ui.js";

document
.getElementById(
"applyCombatBtn"
)
.addEventListener(
"click",
async()=>{

const status =
document.getElementById(
"combatStatus"
);

status.innerHTML =
"Running ComBat...";

try{

const batches =
getBatchAssignments();

const corrected =
await runCombat(

analysisResults.logMatrix,

analysisResults.metadata,

batches

);

analysisResults.correctedMatrix =
corrected;

status.innerHTML =
"ComBat completed";

render();

}
catch(error){

status.innerHTML =
"ComBat failed";

console.error(error);

}

});

let analysisResults;

document
.getElementById("fileInput")
.addEventListener(
"change",
async(e)=>{

const file=e.target.files[0];

analysisResults =
await loadExcel(file);

buildBatchEditor(
analysisResults.metadata
);
  
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
