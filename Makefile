all:sub_yyzq sub_ddjc sub_tztq sub_yysb

SUBDIR1 = ./yyzq
MAKE1 = make
sub_yyzq:
	cd $(SUBDIR1) && $(MAKE1) 

SUBDIR2 = ./ddjc
MAKE2 = make
sub_ddjc:
	cd $(SUBDIR2) && $(MAKE2) 

SUBDIR3 = ./tztq
MAKE3 = make
sub_tztq:
11
	cd $(SUBDIR3) && $(MAKE3) 

SUBDIR4 = ./yysb
MAKE4 = make
sub_yysb:
	cd $(SUBDIR4) && $(MAKE4) 
