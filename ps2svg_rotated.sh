##############################
# ps -> pdf convertion
##############################
ps2pdf $1.ps

##############################
# rotate pdf ccw
##############################
# https://unix.stackexchange.com/questions/394065/command-line-how-do-you-rotate-a-pdf-file-90-degrees
pdftk $1.pdf cat 1-endeast output $1_rotated.pdf

##############################
# pdf -> svg convertion
##############################
pdf2svg $1_rotated.pdf $1_rotated.svg
