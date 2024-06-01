#!/usr/bin/bash
go build -o db_extract *.go
cp db_extract $HOME/bin
