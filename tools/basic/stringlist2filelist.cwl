cwlVersion: v1.0
class: ExpressionTool

label: stringlist2filelist
doc: From a list of file names prefix return a list of files using the extension

requirements:
  - class: InlineJavascriptRequirement
  - class: SchemaDefRequirement
    types:
      - fields:
          - doc: Fastq file for read 1
            name: read_1
            type: File
          - doc: Fastq file for read 2
            name: read_2
            type: File
        name: fastq
        type: record



hints:
  - $import: ubuntu.yml

inputs:
  stringlist:
    type: string[]
  extension:
    type: string
  file_dir:
    type: Directory

outputs:
  files:
    type: File[]

expression:
  "${
      var files = [];
      var l = inputs.file_dir.listing;
      var n = l.length;
      var k = inputs.stringlist.length;
      for (var i = 0; i < n; i++) {
        for (var j = 0; j < k; j++){
          if (l[i].basename.startsWith(inputs.stringlist[j]) && l[i].basename.endsWith(inputs.extension)) {
            files.push(l[i]);
          }
        }
      }
      return { 'files': files};

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

