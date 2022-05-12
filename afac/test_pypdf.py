from PyPDF2 import PdfFileWriter, PdfFileReader
import os, sys

inputpdf = PdfFileReader(open(sys.argv[1], "rb"))
output_pdf, file_extension = os.path.splitext(sys.argv[2])
for i in range(inputpdf.numPages):
    inputpdf = PdfFileReader(open(sys.argv[1], "rb"))
    output = PdfFileWriter()
    output.addPage(inputpdf.getPage(i))
    with open(output_pdf + "-{}".format(i) + file_extension, "wb") as outputStream:
        output.write(outputStream)
