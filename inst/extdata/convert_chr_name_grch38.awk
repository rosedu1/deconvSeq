#!/bin/awk -f
BEGIN{}

{
if ($1~/^MT/)
{gsub(/MT/,"chrM"); print;}
else if ($1~/^ERCC/)
{}
else
{print "chr"$0;}

}

END {}


#awk '{
#if ($1~/^gl/)
#split($1,oldchr,"_");
#$1="";
#print toupper(oldchr[2])".1 "$0;
#}'


