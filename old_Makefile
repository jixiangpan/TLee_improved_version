clean:clean_exe_read_TLee_v07
	-rm -f  *~	

clean_exe_read_TLee_v07:
	-rm exe_read_TLee_v07	
exe_read_TLee_v07:clean_exe_read_TLee_v07
	g++ -o  exe_read_TLee_v07 exe_read_TLee_v07.cc -O -std=c++11 `root-config --libs` -I. -I/home/xji/data0/software/root_build/include -L/home/xji/data0/software/root_build/lib -lMinuit2
