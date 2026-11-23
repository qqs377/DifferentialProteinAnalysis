export function buildBatchEditor(metadata){

const container =
document.getElementById(
"batchEditor"
);

container.innerHTML = "";

const table =
document.createElement("table");

table.className = "batch-table";

table.innerHTML = `

<thead>
<tr>
<th>Sample</th>
<th>Group</th>
<th>Batch</th>
</tr>
</thead>

<tbody>

${metadata.map(sample => `

<tr>
<td>${sample.sample}</td>
<td>${sample.group}</td>

<td>

<select
class="batch-select"
data-sample="${sample.sample}"
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

document
.getElementById("applyCombatBtn")
.disabled = false;

}
