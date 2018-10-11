#duan dian jian ce 
ddjc: c_ddjc.o encode_wav/encode_wav.o
	g++ -g -o $@ $^  -lm

%.o: encode_wav/%.cpp  %.cpp
	gcc -g -c $^

clean:
	rm -rf *.o ddjc  encode_wav/*.o  vad_out_*.wav
