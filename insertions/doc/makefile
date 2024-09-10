TARGET=fcs.pdf

$(TARGET): fcs.xml
	dblatex -P doc.collab.show=0 -P latex.output.revhistory=0 $<

show: $(TARGET)
	evince $<

clean:
	rm $(TARGET)
