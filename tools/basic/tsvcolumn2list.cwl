cwlVersion: v1.0
class: ExpressionTool

label: tsvcolumn2list
doc: Read a TSV table in a file and one column in a list

requirements:
  InlineJavascriptRequirement: {}

hints:
  - $import: ubuntu.yml

inputs:
  table:
    type: File
    inputBinding:
      loadContents: true
  column_name:
    type: string

outputs:
  rows:
    type: string[]

expression:
  "${
     var lines = inputs.table.contents.split('\\n');
     var header = lines[0].split('\\t');
     var colIndex = -1;
     var rows = [];
     for (var i = 0; i < header.length; i++) {
        if (header[i] == inputs.column_name){
           colIndex = i;
           break;
        }
     }
     if (colIndex !== -1){
        for (var i = 1; i < lines.length; i++) {
           var col = lines[i].split('\\t')[colIndex];
           if (col != undefined && col.length != 0){
              rows.push(col)
           }
        }
     }
     return { 'rows': rows } ;
  }"

s:author:
  - class: s:Person
    s:identifier: https://orcid.org/0000-0002-4108-5982
    s:email: mailto:r78v10a07@gmail.com
    s:name: Roberto Vera Alvarez

s:license: https://spdx.org/licenses/OPL-1.0

$namespaces:
  s: http://schema.org/

$schemas:
  - http://schema.org/docs/schema_org_rdfa.html

