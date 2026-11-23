let pyodide;

export async function initializePyodide(){

if(pyodide)
return pyodide;

pyodide =
await loadPyodide();

await pyodide.loadPackage([
"numpy",
"pandas",
"scipy"
]);

return pyodide;

}
