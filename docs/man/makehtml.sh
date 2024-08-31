#!/bin/sh
cp cballs.m cballs.1
man -M . ./cballs.1 | man2html -topm 0 -botm 0 -cgiurl \$title.html > ./cballs.html
