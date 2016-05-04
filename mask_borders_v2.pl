#! /usr/bin/perl

=begin comment

/ DESCRIPTION / ----------------------------------------------------------------

Creates a mask to exclude pixels near the borders of an input image.

First, the borders of the image are traced, and a contour is created using the
ds9 contour function.
The contour converted into a polygon(s), which is subsequently "shrunk" to
exclude the pixel on the borders.
This new polygon(s) is then "filled" with pixels, and finally converted into a
mask.

The user can optionally choose how much to "shrink" the borders and how much to 
"smooth" them (see USAGE).
Polygons with very small numbers of vertices (default = 2 ) can be excluded to
avoid confusion in the "filling" process (see USAGE). 

/ REQUIREMENTS / ---------------------------------------------------------------

> IRAF: IRAF must be active (e.g., source your SciSoft configuration file)

> ds9: if the command "ds9" is not the default in the terminal, setup the
       relevant path in section "CHECK THESE VARIABLES"

> PERL modules:

  Math::Clipper
  Math::Polygon::Tree
  Data::Dumper
  List::MoreUtils
  
  They can be esily insatalled through CPAN, e.g.:
  % su
  % CPAN
  > install Math::Clipper

/ USAGE / ----------------------------------------------------------------------

> Input:

  USAGE:   mask_borders.pl </path/to/file.fits> [-scale <value>] [-sample <value>] [-minv <value>] [-verbose] 
  PRODUCT: borders_<input_filename>

> Parameters:

  o [optional] scale ( default = -10 ):
  
    sets how much the contour encompassing the image bordere are to be scaled
    (it roughly indicates a percentage).
    Following the Math::Clipper convention:
       scale < 0 ==> contour is shrunk
       scale > 0 ==> contour is expanded

  o [optional] sample ( default = 5 ):

    sets how finely the vertices are to be sampled [1 = no sampling].
    Sampling will result in a "smoother" mask.

  o [optional] minv ( default = 2 ):

    sets minimum number of vertices (after rescaling and sampling) for a
    polygon to be used in the mask.
   
  o [optional] verbose ( default = _off_ ):

    reports about execution progress
   
  o [optional] clean ( default = _off_ ):

    removes all the intermediate files (keeping only hte mask)
       
/ EXAMPLE / --------------------------------------------------------------------

> Create a mask for the image file.fits, by shrinking its borders of ~10%, and
  smoothing it by sampling 1 out of 3 pixels:
  
  mask_borders.pl /data/my_project/file.fits -10 3

/ HISTORY / --------------------------------------------------------------------

18/04/2016: mask_borders_v0.pl

22/04/2016: mask_borders_v1.pl

> Update: splitting polygons for fillpoly not to be confused
> Fix: boolean input variables re-commented as "simple" variables

25/04/2016: mask_borders_v2.pl

> Update: not using fillpoly anymore (it was drawing pseudo-random empty
          vertical "stripes"), but instead a PERL library to "fill" the polygons

> Update: added polygon number to region file

/ NOTICE / ---------------------------------------------------------------------

> The contours might contain way too many vertices, so that the routine
  "filling" them with pixels might fail (giving a segmentation fault).
  In that case, try with a larger sample_factor to reduce the number of vertices
  (at the price  of lower resolution).

> Some "holes" within polygons will be ignored !
  (this is a limitation of the Math::Polygon::Tree library)

> "Islands" (i.e. polygons inside polygons inside polygons) will be masked.

/ MEMO / -----------------------------------------------------------------------

> Add option to directly load region file (instead of generating it with ds9).

> Avoid "holes" within polygons to be filled.

/ ACKNOWLEDGMENTS / ------------------------------------------------------------

/ AUTHOR / ---------------------------------------------------------------------

Creator: Paolo Bonfini

Modify at your own risk!

--------------------------------------------------------------------------------

=end comment
=cut

use strict;
use warnings;
use Math::Clipper ':all';
use Math::Polygon::Tree;
use Data::Dumper;
use List::Util qw(min max);
use List::MoreUtils q/natatime/; # to use one-liner iterations (e.g. extracting couples from array)
use File::Basename; # to find path to this script

MAIN:

# CHECK THESE VARIABLES! #######################################################

# user paths -------------------------------------------------------------------
my $bin_ds9 = "ds9";
# path to ds9 executable

# paths ------------------------------------------------------------------------
# NOTE: edit at your own risk!

my $path_program = dirname(__FILE__);
   $path_program =~ s/\.//g; 
# path of this script

# default parameters -----------------------------------------------------------
# NOTE: edit at your own risk!

my $default_scale_factor = -10;
# how much the polygons are to be scaled?

my $default_sample_factor = 5;
# how finely the vertices are to be sampled? [1 = no sampling]

my $default_min_vertex = 2;
# minimum valid vertices for a [rescaled] polygon to be used

################################################################################


################################################################################
  (my $path_to_program = $0)   =~ s/\.\///g;
  # will return "name_of_program"

  my $program_name = (split("/",$path_to_program))[-1];

# parsing input ----------------------------------------------------------------

my $USAGE = "
USAGE: $path_to_program </path/to/file.fits> [-scale <value>] [-sample <value>] [-minv <value>] [-verbose] 
       o file.fits:       image file
       o [scale]:   borders shrinking factor [$default_scale_factor]
       o [sample]:  borders sampling factor  [$default_sample_factor]
       o [minv]:    minimum valid vertices   [$default_min_vertex]
       o [verbose]: report execution messages
       o [clean]:   remove intermediate files
";

my $path_to_file;
# path to input image file

# Variables to check whether relevant keywords have been used:
my $verbose = 0;
my $clean   = 0;
my $scale_factor  = 0;
my $sample_factor = 0;
my $min_vertex = 0;

  if ( grep { $_ eq "-verbose" } @ARGV ){$verbose = 1;}
  if ( grep { $_ eq "-clean"   } @ARGV ){$clean   = 1;}

  if(defined $ARGV[0]){
    $path_to_file = $ARGV[0];
  }else{
    print "ERROR: Input file not provided\n";
    print $USAGE;
    exit;
  }

  for(my $a=0;$a<scalar(@ARGV);$a++){
  # a = argument index
  
    if ($ARGV[$a] =~ /-scale/){
      $scale_factor = 1;
      if(defined $ARGV[$a+1]){
        $scale_factor = $ARGV[$a+1];
        if($verbose){print "$program_name" . ":: scale factor = $scale_factor\n";}
      }else{
        print "ERROR: scale factor not provided\n";
        print $USAGE;
        exit;
      }
    }
  
    if ($ARGV[$a] =~ m/-sample/){
      $sample_factor = 1;
      if(defined $ARGV[$a+1]){
        $sample_factor = $ARGV[$a+1];
        if($verbose){print "$program_name" . ":: sample factor = $sample_factor\n";}
      }else{
        print "ERROR: sample factor not provided\n";
        print $USAGE;
        exit;
      }
    }
  
    if ($ARGV[$a] =~ m/-minv/){
      $min_vertex = 0;
      if(defined $ARGV[$a+1]){
        $sample_factor = $ARGV[$a+1];
        if($verbose){print "$program_name" . ":: minimum valid vertices = $min_vertex\n";}
      }else{
        print "ERROR: sample factor not provided\n";
        print $USAGE;
        exit;
      }
    }
  
  } # END of ARGV
  
  
  if($scale_factor == 0){
    $scale_factor  = $default_scale_factor;
    if($verbose){print "$program_name" . ":: using default scale factor [$scale_factor]\n";}
  }
  
  if($sample_factor == 0){
    $sample_factor  = $default_sample_factor;
    if($verbose){print "$program_name" . ":: using default sample factor [$sample_factor]\n";}
  }
  
  if($min_vertex == 0){
    $min_vertex = $default_min_vertex;
    if($verbose){print "$program_name" . ":: minimum valid vertices [$min_vertex]\n";}
  }
#-------------------------------------------------------------------------------

# setting names ----------------------------------------------------------------

my $file = (split("/",$path_to_file))[-1];
# input file name (striped of path)

my $OUT_BLACKWHITE = "blackwhite.fits";
# name of black &d white image (only 1s and 0s)
# http://www3.sympatico.ca/n.rieck/images/dilbert-03-zeros-and-ones.gif
my $OUT_REG            = "contour.reg";
# contour region file
my $OUT_RESIZED_REG    = "resized_" . $OUT_REG;
# resized contour region file (for debugging purposes)
my $OUT_RESIZED_VERTEX = "resized_vertices.txt";
# resized vertices list
my $OUT_PIXELS          = "pixels.txt";
# pixels within the resized vertices
my $OUT_PIXELS_NESTED = "pixels_nested.txt";
# pixels of polygons nested inside other polygons (will be rejected)
my $OUT_MASK_PIXELS        = "pixels_" . $file;
# output mask
my $OUT_MASK_PIXELS_NESTED = "pixels_nested_" . $file;
# output mask for nested polygons
my $OUT_MASK               = "borders_" . $file;
# final output mask
my $OUT_SPLIT_POLYGONS = "split_polygons";
# file containing the commands to split $OUT_RESIZED_VERTEX into a separate file for each polygon (for debugging purposes)
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
my $NAXIS1 = &hselect($path_to_file,"NAXIS1");
my $NAXIS2 = &hselect($path_to_file,"NAXIS2");
#-------------------------------------------------------------------------------

################################################################################

# converting image to "black & white" ------------------------------------------

# Converting every pixel different from 0 to 1:
if(-e $OUT_BLACKWHITE){system("rm $OUT_BLACKWHITE")}
&imexpr("$OUT_BLACKWHITE","( (a!=0) ? 1 : 0 )","$path_to_file");

#-------------------------------------------------------------------------------

# creating contour region file -------------------------------------------------

system("$bin_ds9 $OUT_BLACKWHITE -contour yes -contour nlevels 1 -contour smooth 1 -contour close -contour convert -regions system image -regions save $OUT_REG -exit");

#-------------------------------------------------------------------------------

# loading region file ----------------------------------------------------------

my $units;
# units of input region file [image,physical]

my @all_polygons;
# array of polygons: contains references to arrays containing the refernces to
# the vertices for all of the polygons read from the region file <array>
# 
# $all_polygons = [
#   [
#    [446,1383],[445,1382]
#   ]
#   [
#    [895,695],[894,694],[893,1023]
#   ]
#   ...
#   ]
# ]

my $counter = 0;
# counter for all_plygons


  open(OUT_REG,"<$OUT_REG") or die " ERROR: cannot open region file $OUT_REG";  
  while ($_=<OUT_REG>) {		  	        
    
    if ($_ =~ /image.*/)    {$units = "image"}
    if ($_ =~ /physical.*/) {$units = "physical"}
    
    
    if ($_ =~ /polygon.*/) {
    
      $counter++;
    
      my $line = $_;
         $line =~ s/polygon\(//g;
         $line =~ s/\)//g;
      chomp($line); 
       
      my @read = split(",",$line); 
      my $xy = natatime( 2 , @read );
      # splitting vertices for current polygon into couples

      my $polygon;
      # curent polygon
      #
      # $polygon = [
      #    [446,1383],[445,1382]
      # ]

      # Iterating over vertices coordinates:
      while (my ($x,$y) = $xy->())
      {
	
	# creating reference to vertex:
	my $ref_vertex = [$x,$y];
	
	# creating polygon:
	push(@{$polygon},$ref_vertex);
      }

      push(@all_polygons,$polygon);

    }
   
  }
  close(OUT_REG);
  
#-------------------------------------------------------------------------------

# resizing polygon -------------------------------------------------------------

my $offset_polygons = offset(\@all_polygons, $scale_factor);
# Clipper.offset requires a reference to an array of polygons

# Debugging:
# print Dumper(@all_polygons);
# print Dumper($offset_polygons);
#-------------------------------------------------------------------------------

# creating resized vertex/region file ------------------------------------------
#
# NOTE: the output vertices of the resized polygons will be sampled 1 every
#       $sample_factor

  my @polygon_n_sampled_vertices;
  # number of (sampled) vertices contained in each polygon

  open(OUT_RESIZED_VERTEX,">$OUT_RESIZED_VERTEX");
  open(OUT_RESIZED_REG   ,">$OUT_RESIZED_REG"   );

  print OUT_RESIZED_REG "global color=red dashlist=8 3 width=1 font=\"helvetica 10 normal roman\" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1\n";
  print OUT_RESIZED_REG "$units\n";


  my $p = 0;
  # polygon counter
  
  foreach my $polygon (@{$offset_polygons}){
  # curent polygon
  #
  # $polygon = [
  #    [446,1383],[445,1382]
  # ]

    $polygon_n_sampled_vertices[$p] = 0;

    print OUT_RESIZED_REG "polygon(";
    print OUT_RESIZED_VERTEX "# polygon $p\n";


    for(my $v=0;$v<scalar(@{$polygon});$v+=$sample_factor){
    # v = vertex counter
 
      $polygon_n_sampled_vertices[$p]++;

      # vertex coordinates:
      my $x = $polygon->[$v][0];
      my $y = $polygon->[$v][1];
      
      printf OUT_RESIZED_VERTEX "%-6.0f  %-6.0f\n", $x, $y;
      print OUT_RESIZED_REG    "$x,$y";
    
      if( ($v+$sample_factor)<scalar(@{$polygon}) ){
        print OUT_RESIZED_REG ",";
      }
    }

    print OUT_RESIZED_REG ") # text = \"p $p\"\n";
 
    $p++;
 
 } #END of polygons
  
  close(OUT_RESIZED_VERTEX);
  close(OUT_RESIZED_REG);

#-------------------------------------------------------------------------------

# summary ----------------------------------------------------------------------
  if($verbose){
    printf "$program_name" .":: %s polygon(s) created\n", scalar(@polygon_n_sampled_vertices);
  }
#-------------------------------------------------------------------------------

# filling resized contour with pixles ------------------------------------------

  if($verbose){
    printf "$program_name" .":: filling polygon(s), this will take a while ...\n";
  }

  if(-e $OUT_PIXELS){system("rm $OUT_PIXELS");}
  # IMPORTANT: appending the "filled" pixels for each polygon


  # Splitting polygons, and filling each one separately:
  
  my $line_head = 0;
  my $line_tail = 0;
  # tmp

  my $polygons_n_actual = 0;
  # polygons surviving after screening for min vertices
  my $polygons_n_nested = 0;
  # nested polygons

  open(OUT_SPLIT_POLYGONS,">$OUT_SPLIT_POLYGONS");

  open(OUT_PIXELS       ,">$OUT_PIXELS");
  open(OUT_PIXELS_NESTED,">$OUT_PIXELS_NESTED");

  POLYGONS: for(my $p=0;$p<scalar(@polygon_n_sampled_vertices);$p++){
  # p = polygon counter

    if($verbose){
      my $tot = scalar(@polygon_n_sampled_vertices);
      my $polygon_n = $p + 1;
      system("echo -n \"\\rpolygon [ID=$p] $polygon_n/$tot\"");
    }

    # defining current polygon -------------------------------------------------

    my $nested = 0;
    # is polygon nested inside an other poligon?

    my $polygon = @{$offset_polygons}[$p];
    # NOTE: @polygon_n_sampled_vertices and @{$offset_polygons} should be indexed in the same way.
    #       However it doesn't really matter, since the important is to through all of the polygons

    my $boundary = Math::Polygon::Tree->new( $polygon );
    # polygon boundary object
  
    # rejections ---------------------------------------------------------------
    # Pixels belonging to these polygons will not be filled

    # > Small polygons

    if($polygon_n_sampled_vertices[$p] < $min_vertex){
    
      $line_head += $polygon_n_sampled_vertices[$p] + 1;
      # +1 goes for the "#" separating each vertices group
      $line_tail =  $polygon_n_sampled_vertices[$p];
      next POLYGONS;
    
    }

    # > Polygon inside an other polygon (note: "islands" will be excluded too)

    my $counter_other = 0;
    # tmp: other polygon counter
    
    OTHER_POLYGONS: foreach my $other_polygon (@{$offset_polygons}){
    # every other polygon ...

      $counter_other++;
      
      my $other_boundary = Math::Polygon::Tree->new( $other_polygon );

      my $contains_polygon_rough = $other_boundary->contains_polygon_rough( $polygon );
      
      if ( (defined($contains_polygon_rough)) && ($contains_polygon_rough == 1) )  {
      # conservative: some "undefined" results may actually indicate nested polygons
      
        # Debug:
        #if($verbose){
        #  print "$program_name" .":: polygon $p is inside an other polygon\n";
        #}
  
        $nested = 1;
	$polygons_n_nested++;
	last OTHER_POLYGONS;
      }
    
    }
    

    # splitting polygons -------------------------------------------------------

    $polygons_n_actual++;
  
    my $OUT_RESIZED_VERTEX_POLYGON = "polygon_$p.txt";
    # temp: vertices for polygon $p

    $line_head += $polygon_n_sampled_vertices[$p] + 1;
    # +1 goes for the "#" separating each vertices group
    $line_tail =  $polygon_n_sampled_vertices[$p];

    my $DO_split_polygon = "head -$line_head $OUT_RESIZED_VERTEX | tail -$line_tail > $OUT_RESIZED_VERTEX_POLYGON";

    `$DO_split_polygon`;

    print OUT_SPLIT_POLYGONS "$DO_split_polygon\n";
    # storing commands for debugging
    

    # finding min and MAX coordinates of a polygon -----------------------------
    
    # min and MAX coordinates expected in a polygon (to speed up iterations through pixels):
    my $polygon_x_min = 0;
    my $polygon_x_MAX = $NAXIS1;
    my $polygon_y_min = 0;
    my $polygon_y_MAX = $NAXIS2;
    # [pixel]

    # Polygon vertices coordinates:
    my @polygon_vertices_x;
    my @polygon_vertices_y;
    # [pixel]

    open(OUT_RESIZED_VERTEX_POLYGON,"<$OUT_RESIZED_VERTEX_POLYGON") or die " ERROR: cannot open polygon file $OUT_RESIZED_VERTEX_POLYGON";  
    while ($_=<OUT_RESIZED_VERTEX_POLYGON>) {		  	        
    
      if($_ !~ /^#.*/){
         my @read = split(" ",$_);
      
         push(@polygon_vertices_x,$read[0]);
         push(@polygon_vertices_y,$read[1]);
      }
   
    }
    close(OUT_RESIZED_VERTEX_POLYGON);
    
    $polygon_x_min = min(@polygon_vertices_x);
    $polygon_x_MAX = max(@polygon_vertices_x);
    $polygon_y_min = min(@polygon_vertices_y);
    $polygon_y_MAX = max(@polygon_vertices_y);
    
    # filling polygon with pixels ----------------------------------------------

    for(my $x=$polygon_x_min;$x<$polygon_x_MAX;$x++){
    # x = x coordinate index
    
      for(my $y=$polygon_y_min;$y<$polygon_y_MAX;$y++){
      # y = y coordinate index
 
        my $point = [$x,$y];

        my $is_inside = $boundary->contains( $point );

        if ( $is_inside ) {
	
          if($nested == 0) {printf OUT_PIXELS        "%-10s %-10s\n", @{$point}[0], @{$point}[1];}
          if($nested == 1) {printf OUT_PIXELS_NESTED "%-10s %-10s\n", @{$point}[0], @{$point}[1];}
        }	
     
      }
    
    }

      
    system("rm $OUT_RESIZED_VERTEX_POLYGON");
  }

  close(OUT_SPLIT_POLYGONS);
  close(OUT_PIXELS);
  close(OUT_PIXELS_NESTED);
  
  if($verbose){
    system("echo \"\"");
  }
  
#-------------------------------------------------------------------------------

# creating mask for valid pixels -----------------------------------------------
  if($verbose){
    print "$program_name" .":: generating mask, be patient ...\n";
  }

  if(-e $OUT_MASK_PIXELS){system("rm $OUT_MASK_PIXELS");}
  &badpiximage("$OUT_PIXELS","$path_to_file",$OUT_MASK_PIXELS,1,0);
#-------------------------------------------------------------------------------

# creating mask for nested pixels ----------------------------------------------
  if($verbose){
    print "$program_name" .":: generating mask for nested pixels, be patient ...\n";
  }

  if(-e $OUT_MASK_PIXELS_NESTED){system("rm $OUT_MASK_PIXELS_NESTED");}
  &badpiximage("$OUT_PIXELS_NESTED","$path_to_file",$OUT_MASK_PIXELS_NESTED,0,1);
#-------------------------------------------------------------------------------

# combining masks --------------------------------------------------------------
if(-e $OUT_MASK){system("rm $OUT_MASK")}
&imexpr("$OUT_MASK","( (a && (!b)) ? 1 : 0 )","$OUT_MASK_PIXELS","$OUT_MASK_PIXELS_NESTED");
#-------------------------------------------------------------------------------


# summary ----------------------------------------------------------------------
  if($verbose){
    printf "$program_name" .":: %-4s [/%s] masked   polygon(s) \n", ($polygons_n_actual - $polygons_n_nested), scalar(@polygon_n_sampled_vertices);
    printf "$program_name" .":: %-4s [/%s] nested   polygon(s) \n", $polygons_n_nested, scalar(@polygon_n_sampled_vertices);
    printf "$program_name" .":: %-4s [/%s] rejected polygon(s) (< %s vertices)\n", $polygons_n_nested, scalar(@polygon_n_sampled_vertices), $min_vertex;
  }
#-------------------------------------------------------------------------------

# cleaning ---------------------------------------------------------------------
  if($clean == 1){
    if($verbose){
      print "$program_name" .":: removing intermediate files\n";
    }
  
    system("rm $OUT_BLACKWHITE");
    system("rm $OUT_REG");
    system("rm $OUT_RESIZED_REG");
    system("rm $OUT_RESIZED_VERTEX");
    system("rm $OUT_SPLIT_POLYGONS");
    system("rm $OUT_PIXELS");
    system("rm $OUT_PIXELS_NESTED");
    system("rm $OUT_MASK_PIXELS");
    system("rm $OUT_MASK_PIXELS_NESTED");

  } # END of cleaning
#-------------------------------------------------------------------------------
################################################################################
# EOF


################################################################################
# FUNCTIONS ####################################################################
################################################################################

# hselect driver ---------------------------------------------------------------
# USAGE: hselect(path_to_file,keyword)
# RETURN: [keyword_value,"ERR"]

sub hselect{

  my $keyword_value;

  my $path_to_file = $_[0];
  # input file name

  my $filename = (split("/",$path_to_file))[-1];
  # input file name (striped of path)

  open(HSELECT,">hselect_$filename.cl");
  print HSELECT <<EOF;
images
imutil
unlearn hselect
hselect images=$path_to_file fields="\$I,$_[1]" expr=yes
logout
EOF

  close(HSELECT);

  system("cl < hselect_$filename.cl > hselect_$filename.log");

  if(system("grep \"$filename\" hselect_$filename.log >/dev/null") == 0){
 
    my $DO_awk = "awk '{if (\$0 ~ /\Q$filename\E/) printf \"\%s\\n\", \$4}' hselect_$filename.log";
    # capturing the 4th field of line containing file name
    # NOTE: this depends on the hselect output

    $keyword_value = `$DO_awk`;
    chomp $keyword_value;
  }

  system("rm hselect_$filename.cl hselect_$filename.log");
  if(-e "uparmimlhselet.par"){
    system("rm uparmimlhselet.par");
  }

  if ($keyword_value){
    return $keyword_value;
  }else{
    return "ERR";
  }

}
#-------------------------------------------------------------------------------

# imexpr driver  ---------------------------------------------------------------
# USAGE: imexpr(output.fits,"expression",a.fits,b.fits,..,z.fits)

sub imexpr{

  my $n_operators = scalar(@_) - 2;
  # number of operators

  my @operator_labels = ("a".."z");
  
  my $operators_string = "";
  
  for(my $i=2;$i<scalar(@_);$i++){
  # i = input index
    $operators_string .= " " . "$operator_labels[$i-2]=" . "$_[$i]";
  }

  open(IMEXPR,">imexpr_$_[0].cl");
  print IMEXPR <<EOF;
images
imutil
unlearn imexpr
imexpr expr="$_[1]" output=$_[0] $operators_string
logout
EOF

  close(IMEXPR);

  system("cl < imexpr_$_[0].cl > imexpr_$_[0].log");
  system("rm imexpr_$_[0].cl imexpr_$_[0].log");
  if(-e "uparmimlimexpr.par"){
    system("rm uparmimlimexpr.par");
  }
}
#-------------------------------------------------------------------------------

# badpixim driver  -------------------------------------------------------------
# USAGE: badpiximage(pixels.txt,template.fits,out.fits,[goodvalue],[badvalue])

sub badpiximage{

  my $goodvalue = 1;
  my $badvalue  = 0;
  
  my $n_parameters = scalar(@_);
  # number of input parameters
  
  if($n_parameters > 3){

    if ( ($n_parameters == 4) || ($n_parameters > 5) ) {
      print "badpixim:: expecting 3 or 5 parameters";
      return "ERR";
    }else{
     $goodvalue = $_[3];
     $badvalue  = $_[4];
    }
  }

  open(BADPIXIMAGE,">badpiximage_$_[0].cl");
  print BADPIXIMAGE <<EOF;
noao
imred
ccdred
unlearn badpiximage
badpiximage fixfile="$_[0]" template="$_[1]" image="$_[2]" goodvalue=$goodvalue badvalue=$badvalue
logout
EOF

  close(BADPIXIMAGE);

  system("cl < badpiximage_$_[0].cl > badpiximage_$_[0].log");
  system("rm badpiximage_$_[0].cl badpiximage_$_[0].log");
  if(-e "uparmccdbadpie.par"){
    system("rm uparmccdbadpie.par");
  }
}
#-------------------------------------------------------------------------------
################################################################################

