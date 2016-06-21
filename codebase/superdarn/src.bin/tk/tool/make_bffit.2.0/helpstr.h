/* hlpstr.h
   ========
   Author: R.J.Barnes
*/

/*
 $License$
*/


char *hlpstr[]={
"make_bffit - Creates a bffit fitacf format file from a rawacf format file.\n",
"make_bffit --help\n",
"make_bffit rawname fitname\n",
"make_bffit -new [rawacfname] > fitname\n",

"--help\tprint the help message and exit.\n",
"rawname\tfilename of the raw (dat) format file.\n",
"fitname\tfilename of the fit format file to create.\n",
"-new\tinput file is in  rawacf file format and the output is in fitacf file format.\n",
"rawacfname\tfilename of the rawacf format file. If this is omitted the file is read from standard input.\n",

NULL};
