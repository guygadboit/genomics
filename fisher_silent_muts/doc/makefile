TARGET=silent_muts.pdf

$(TARGET): silent_muts.xml
	dblatex -P doc.collab.show=0 -P latex.output.revhistory=0 $<

show: $(TARGET)
	evince $<

clean:
	rm $(TARGET)
