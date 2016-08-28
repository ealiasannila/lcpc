//-------------------------------------------------------------------------/
//Copyright (C) 2003, 2004, 2005, ALL RIGHTS RESERVED.
//Centre for Sys. Eng. & App. Mech.           FEMAGSoft S.A.
//Universite Cathalique de Louvain            4, Avenue Albert Einstein
//Batiment Euler, Avenue Georges Lemaitre, 4  B-1348 Louvain-la-Neuve
//B-1348, Louvain-la-Neuve                    Belgium
//Belgium
//-------------------------------------------------------------------------/
//
//Name:         parse.h ( parse command lines options with help and copyright
//                        information )
//Author:       Liang, Wu (wu@mema.ucl.ac.be, wuliang@femagsoft.com)
//Created:      11/2005
//-------------------------------------------------------------------------/

#ifndef  _PARSE_H
#define  _PARSE_H

#include <string>
#include <iostream>
#include <stdio.h>      /* printf, fopen */
#include <stdlib.h>     /* exit, EXIT_FAILURE */


class Options
{
   public:
     Options():parsebdm(true),showme(false),tecplot(false),metapost(true), debug(false), fileindex(1) 
     { 
     }
     
     bool parsebdm, showme, tecplot, metapost, debug;
     int  fileindex;
};

void help()
{
	std::cout<<"Poly2Tri: Fast polygon triangulation by sweep line algorithm.\n";
    std::cout<<"Copyright (C) 2003, 2004, 2005, ALL RIGHTS RESERVED.\n";
    std::cout<<"Centre for Sys. Eng. & App. Mech.           FEMAGSoft S.A.\n";
    std::cout<<"Universite Cathalique de Louvain            4, Avenue Albert Einstein\n";
    std::cout<<"Batiment Euler, Avenue Georges Lemaitre, 4  B-1348 Louvain-la-Neuve\n";
    std::cout<<"B-1348, Louvain-la-Neuve                    Belgium\n";
    std::cout<<"Belgium\n\n";
    std::cout<<"Usage: poly2tri [-pstmdh] filename.\n";
    std::cout<<"Options:\n";
    std::cout<<"     -p: DO NOT parse input bdm (boundary mesh) file.\n";
    std::cout<<"         Use this option if the input is very large without any comments.\n";
    std::cout<<"     -s: Output results as ShowMe formats.\n";
    std::cout<<"     -t: Output results as TECPLOT ASCII PLT format.\n";
    std::cout<<"     -m: DO NOT output monotone polygons and triangulation results\n";
    std::cout<<"         as metapost source file.\n";
    std::cout<<"     -d: Output debug information (log file).\n";
    std::cout<<"     -h: This help message.\n\n";
    std::cout<<"Examples:\n";
    std::cout<<"poly2tri sample.bdm\n";
    std::cout<<"     parse and triangulate sample.bdm and output result as metapost source file.\n";
    std::cout<<"poly2tri -pm sample.bdm\n";
    std::cout<<"     triangulate sample.bdm and no result output!.\n";
    std::cout<<"poly2tri crzaybox.bdm -T -D\n";
    std::cout<<"     triangulate crazybox.bdm and no result output!\n";
    std::cout<<"Report bugs to wu@mema.ucl.ac.be or wuliang@femagsoft.com.\n\n";
}

void parse(int argc, char **argv, Options& options)
{
  unsigned int num=5;
  unsigned int files=0, fileindex=1;
  
  std::string  option, start, end;
  
  if( argc<2 ) { help(); exit(1); }
  for (int k = 1; k < argc; k++)
  {
	 option=argv[k];
	 if(option[0]!='-') { fileindex=k; files++; continue; } 
	 if(option.size() >= num+1 || option[0]!='-' || option.size() <= 1 ) { help(); exit(1); }
         for(unsigned int i=1; i < option.size(); i++) 
	 {
		 switch (option[i])
		 {
			 case 'P': options.parsebdm=false; break;
			 case 'p': options.parsebdm=false; break;
			 case 'S': options.showme=true;    break;
			 case 's': options.showme=true;    break;
			 case 'T': options.tecplot=true;   break;
			 case 't': options.tecplot=true;   break;
			 case 'M': options.metapost=false; break;
			 case 'm': options.metapost=false; break;
			 case 'D': options.debug=true;     break;
			 case 'd': options.debug=true;     break;
		         case 'h': help(); exit(1);        break;
			 case 'H': help(); exit(1);        break;
			 default : help(); exit(1);        break;
		 }
	 }
   }

  options.fileindex=fileindex;
  
  if( files > 1 ) { std::cout<<"only one file, please!\n"; exit(1); }
  else if (files <=0 ) { help(); exit(1); }

}

#endif

