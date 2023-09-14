const form = document.getElementById("form");
const fileInput = form.elements['sdf'];
const submitButton = document.getElementById("submit");
const output = document.querySelector("output");
const outputTableBody = document.getElementById("output-table-body");
const spinner = document.getElementById("spinner");
const drawnMolIgnoredMsg = document.getElementById("drawn-mol-ignored-msg");
const error = document.getElementById("error");

let drawnMoleculePresent = false;
let inputFilePresent = false;

async function submitMolecules()
{
  if (!inputFilePresent && !drawnMoleculePresent)
    return;

  setOutputVisible(false);
  setError(null);

  const formData = new FormData(form);

  if (drawnMoleculePresent && !inputFilePresent)
    formData.set("sdf", new Blob([JSDraw.get("jsdraw").getMolfile()]));

  try
  {
    setSpinnerVisible(true);

    const res = await
      fetch("/deep-dili/predict-sdf", {
        method: "POST",
        body: formData,
        headers: { 'Accept': 'application/json' }
      });

    if (res.ok)
      setResults(await res.json());
    else
      setError(await res.json());
  }
  finally
  {
    setSpinnerVisible(false);
  }
}

function setResults(results)
{
  while(outputTableBody.firstChild)
    outputTableBody.firstChild.remove();

  for (const result of results)
    outputTableBody.appendChild(createResultTableRow(result));

  setOutputVisible(true);
}

function createResultTableRow(result)
{
  const tr = document.createElement("tr");

  const molIdTd = document.createElement("td");
  molIdTd.textContent = result.mol_id;
  tr.appendChild(molIdTd);

  const predTd = document.createElement("td");
  predTd.textContent = result.prob_pred.toFixed(2);
  tr.appendChild(predTd);

  const classTd = document.createElement("td");
  const toxic = result.class_pred >= 0.5;
  classTd.textContent = toxic ? "DILI" : "non-DILI";
  classTd.classList.add(toxic ? "toxic" : "non-toxic");
  tr.appendChild(classTd);

  return tr;
}

function setOutputVisible(visible)
{
  output.style.display = visible ? "block" : "none";
}

function setError(errorMessage)
{
  if (errorMessage)
  {
    if (typeof errorMessage == 'string')
      error.textContent = errorMessage;
    else
      error.textContent = "The input could not be processed. Please check the format of the inputs.";

    error.style.display = "block";
  }
  else
    error.style.display = "none";
}

function setDrawnMoleculeIgnoredWarningVisible(visible)
{
  drawnMolIgnoredMsg.style.display = visible ? "inline" : "none";
}

function setSpinnerVisible(visible)
{
  spinner.style.display = visible ? "block" : "none";
}

function fileInputChanged()
{
  inputFilePresent = fileInput.value.length !== 0;
  inputChanged();
}

function drawnMoleculeChanged(jsdraw)
{
  drawnMoleculePresent = jsdraw.getMolWeight() > 0;
  inputChanged();
}

function inputChanged()
{
  submitButton.disabled = !drawnMoleculePresent && !inputFilePresent;
  setOutputVisible(false);
  setError(null);
  setDrawnMoleculeIgnoredWarningVisible(drawnMoleculePresent && inputFilePresent);
}

fileInput.addEventListener('change', fileInputChanged);

submitButton.addEventListener('click', function() {
  if (drawnMoleculePresent || inputFilePresent)
    submitMolecules();
});

inputChanged();