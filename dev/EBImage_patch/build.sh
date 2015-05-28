#!/bin/bash

R CMD Sweave changes.Rnw && pdflatex changes.tex && pdflatex changes.tex
