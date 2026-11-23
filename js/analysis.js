export async function loadExcel(file){

const data =
await file.arrayBuffer();

const workbook =
XLSX.read(data);

const sheet =
workbook.Sheets[
workbook.SheetNames[0]
];

const json =
XLSX.utils.sheet_to_json(sheet);

const abundanceCols = Object.keys(json[0])
.filter(c =>
c.includes(
"Abundances (Normalized)"
)
);

const sampleInfo =
abundanceCols.map(col=>{

const match =
col.match(
/:\s*([^:]+):\s*Sample,\s*(.+)$/
);

return{

column:col,
sample:match[1],
group:match[2]

};

});

return {

raw:json,

metadata:sampleInfo,

pca:null,

volcano:null,

heatmap:null

};

}
