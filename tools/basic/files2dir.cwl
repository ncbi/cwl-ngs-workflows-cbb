cwlVersion: v1.0
class: ExpressionTool

label: files2dir
doc: Group all input files in a directory

requirements:
  InlineJavascriptRequirement: {}

hints:
  - $import: ubuntu.yml

inputs:
  files:
    type: File[]
  dir:
    type: string

outputs:
  output:
    type: Directory

expression: |
  ${
     var listing = [];
     for (var i = 0; i < inputs.files.length; i++) {
        listing.push(inputs.files[i]);
     }
     return {
        "output": {
          "class": "Directory",
          "basename": inputs.dir,
          "listing": listing
        }
     };
  }
