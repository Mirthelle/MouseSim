#!/usr/bin/perl -w 

use strict; 
use LWP::Simple;

use IO::File;

use GO::TermFinder;
use GO::AnnotationProvider::AnnotationParser;
use GO::OntologyProvider::OboParser;

use GO::Utils::File    qw (GenesFromFile);
use GO::Utils::General qw (CategorizeGenes);

###################################################################################
sub Usage{
###################################################################################

    my $message = shift;

    if (defined $message){

		print $message, "\n";
    }

    print <<USAGE;

##############################################
##                                          ##
##              M o u s e S i m             ##
##                                          ##
##############################################


This program takes a PFAM_ID as input and will fetch a FASTA file
with human sequences part of the selected family. Then it will get
GO terms associated with the list of genes and find similar genes 
annotated to the same GO nodes in mouse genome. 

In addition to PFAM_ID you also can provide names for 3 extra 
files: a human annotation file (with GO terms), a mouse annotation 
file (with GO terms) and a .obo file wich will contain the 
ontology. If you only provide a PFAM_ID, default annotation and 
ontology files will be used (you can find them in the subdirectory
"files").

This will generate an output files with the results named "mouseSim_results.txt".

Usage:

\t ./mouseSim.pl <PFAM_ID> <human_annotation_file> <mouse_annotation_file> <obofile> 

e.g.

\t ./mouseSim.pl PF00870 ../t/gene_association.sgd ../gene_association.mgi ../t/gene_ontology_edit.obo

USAGE

    exit;

}


######################################################################
sub getFASTA {                                                      ##
######################################################################
## Gets FASTA from PFAM url using PFAM accession ID as input.       ##
######################################################################
    my($acc) = @_;
    
    # Getting fasta entry
    my $loc="http://pfam.xfam.org/family/"; 
    my $data = "/alignment/seed/format?format=fasta&alnType=seed&order=a&case=l&gaps=none&download=0";
    my $all = $loc . $acc . $data;
    my $entry = get($all);
    my @fasta = split("\n", $entry);    
}

######################################################################
sub humanGenesNames {                                               ##
######################################################################
## Using a fasta file as input writes a text file with human        ##
## gene names.                                                      ## 
######################################################################
    my(@fasta) = @_;
    
    print "-- Selecting human gene names from FASTA\n";
   
    # Selecting headers correspondig to human genes and storing them in @headers.
    my @headers = ();
    foreach my $line (@fasta) {
        if (($line =~  m/^>/)  and ($line =~ m/_HUMAN/)) {
            my $header = $line;
            push @headers, $header;
        }
    }
    
    # Selecting only gene names and storing them in @gene_names.
    my @gene_names = ();
    foreach my $name (@headers) {
        $name =~ s/>//;     # Remove first >
        my @parts = split("_", $name);  # Split header in "_"
        push @gene_names, $parts[0];    # Select only first element which is gene name
    }
    print "-- Human genes selected were: ", "@gene_names", "\n";


    # Writing text file with one gene name per line.
    my $filename = "names.txt";
    open(my $fh, '>', $filename) or die "Could not open file '$filename' $!";
    print $fh join("\n", @gene_names);
    close $fh;
    return $filename, [@gene_names];
}

######################################################################
sub GOterms {                                                       ##
######################################################################
## Takes 3 file names as inputs: gene_list, annotation_file         ##
## (annotated genome) and obo_file (ontology_file).                 ##
## Returns hash with key being aspects (P, C, F) and values being   ##
## arrays with annotated GO terms in that aspect.                   ##
######################################################################
    my($annotation_file, $obo_file, $genes_list) = @_;
    
    my $total_num = 46661; # Number of genes in annotated genome.
    
    ## Setting up the objects we need
    my $process   = GO::OntologyProvider::OboParser->new(ontologyFile => $obo_file,
                                 aspect       => 'P');
    my $component = GO::OntologyProvider::OboParser->new(ontologyFile => $obo_file,
                                 aspect       => 'C');
    my $function  = GO::OntologyProvider::OboParser->new(ontologyFile => $obo_file,
                                 aspect       => 'F');
    
    my $annotation = GO::AnnotationProvider::AnnotationParser->new(annotationFile=>$annotation_file);
    
    my $termFinderP = GO::TermFinder->new(annotationProvider=> $annotation,
                          ontologyProvider  => $process,
                          totalNumGenes     => $total_num,
                          aspect            => 'P');
    
    my $termFinderC = GO::TermFinder->new(annotationProvider=> $annotation,
                          ontologyProvider  => $component,
                          totalNumGenes     => $total_num,
                          aspect            => 'C');
    
    my $termFinderF = GO::TermFinder->new(annotationProvider=> $annotation,
                          ontologyProvider  => $function,
                          totalNumGenes     => $total_num,
                          aspect            => 'F');
        
    my $cutoff = 0.05;
    
    ## Gene categorization (and error control)
        my @genes = GenesFromFile($genes_list);
        unlink $genes_list;
        my (@list, @notFound, @ambiguous);
    
        CategorizeGenes(annotation  => $annotation,
                genes       => \@genes,
                ambiguous   => \@ambiguous,
                unambiguous => \@list,
                notFound    => \@notFound);
    
         if (@list){
             print "-- Fetching standard DB gene names: \n";
             foreach my $gene (@list){
                 print "\t", $gene, "\t", $annotation->standardNameByName($gene), "\n";
             }
    
             print "\n";
         }
         else{
             print "None of the gene names were recognized\n";
             print "They were:\n\n";
             print join("\n", @notFound), "\n";
             next;
         }
    
         if (@ambiguous){
             print "The following gene(s) are ambiguously named, and so will not be used:\n";
             print join("\n", @ambiguous), "\n\n";
         }
    
         if (@notFound){
             print "The following gene(s) were not recognized, and will not be considered:\n\n";
             print join("\n", @notFound), "\n\n";
         }
        
        ## Getting GO terms for each aspect and storing them into an array of arrays.
        my @aspects;
        foreach my $termFinder ($termFinderP, $termFinderC, $termFinderF){  
            my @pvalues = $termFinder->findTerms(genes => \@list, calculateFDR => 1);
            
            my @goterms = ();
            for (my $i=0; $i<(scalar @pvalues); $i++) {
                my $goid = $pvalues[$i]{"NODE"}{"GO::Node::__goid"};
                push @goterms, $goid;       
            }
            push @aspects, [@goterms];
         }
        
        ## Storing each array in aspects into a hash, key being the aspect.
         my %GOterms = (
            'P' => [@{$aspects[0]}],
            'C' => [@{$aspects[1]}],
            'F' => [@{$aspects[2]}]
         );
         
         return %GOterms;
}

######################################################################
sub findSimilarGenes {
######################################################################
## Given the name of a file with annotated genome for mouse and a   ##
## hash associating aspects with GO IDs, returns mouse genes which  ##
## are anotated to the same GO ID nodes given as input. The result  ##
## is a hash of hashes, first key being the aspect and the second   ##
## one a GO ID which has as value an array with gene names.         ##
######################################################################
    my($mouse_anno_file, %human_terms) = @_;
    
    open(ANNOFILE, $mouse_anno_file);
    
    ## Saving in %mouse_terms GO terms as key and arrays with lists of genes as value.
    my %mouse_terms;
    
    foreach my $aspect (keys %human_terms){
        
        my %mouse_GO;
        foreach my $line (<ANNOFILE>) {
            if ($line =~ m/\t$aspect\t/) {
                my @GO_line = split (" ", $line);
                my @gene_line = split (" ", $line);
                
                if (exists $mouse_GO{$GO_line[3]}) {
                    push $mouse_GO{$GO_line[3]}, $gene_line[2]
                }
                else {
                    $mouse_GO{$GO_line[3]} = [];
                    push $mouse_GO{$GO_line[3]}, $gene_line[2];
                }
            }
        }
        $mouse_terms{$aspect} = {%mouse_GO};
    }
    close ANNOFILE;
 
    ## Finding genes in mouse with the same GO ID than the input.       
    my %similar_genes;
    foreach my $aspect (keys %human_terms){
        my %GOandGenes;
        for (my $i=0; $i < (scalar @{$human_terms{$aspect}}); $i++) {
            my $GO_term = $human_terms{$aspect}[$i];
            if (exists $mouse_terms{$aspect}{$GO_term}) {
                $GOandGenes{$GO_term} = $mouse_terms{$aspect}{$GO_term};
            }            
        }
        $similar_genes{$aspect} = {%GOandGenes};
    }
    return %similar_genes;
}

######################################################################
sub moreSimilarGenes {                                              
######################################################################
## Using the output of findSimilarGenes, counts how many times each ##
## gene appears to obtain which ones are more similar to human.     ##
######################################################################
    my(%similar_genes) = @_;
    
    ## Getting frequency for each gene
    my %genes_freq;
    foreach my $aspect (sort keys %similar_genes) {
        foreach my $GO_ID (sort keys %{$similar_genes{$aspect}}) {
            for (my $i = 0; $i < (scalar @{$similar_genes{$aspect}{$GO_ID}}); $i++) {
                my $gene = $similar_genes{$aspect}{$GO_ID}[$i];
                if (exists $genes_freq{$gene}) {
                    $genes_freq{$gene}++; 
                }
                else{
                    $genes_freq{$gene} = 1;
                }
            }
        }
    }
    return %genes_freq;
}

##############################################
### 	   M A I N 	  P R O G R A M         ##
##############################################
{
	&Usage if (@ARGV < 1);

	my ($PFAM_ID, $anno_human_file, $anno_mouse_file, $obo_file);

	if (@ARGV == 1) {
		$PFAM_ID = $ARGV[0];
		$anno_human_file = "files/gene_association.goa_human";
		$anno_mouse_file = "files/gene_association_mouse.mgi";
		$obo_file = "files/go-basic.obo";

		print "-- No annotation files nor ontology file specified, using default. \n"; 
	}

	if (@ARGV == 4) {
		$PFAM_ID = shift;
		$anno_human_file = shift;
		$anno_mouse_file = shift;
		$obo_file = shift;
	}  

    if ($PFAM_ID =~ m/^PF/) {
        print "-- Starting analysis\n";

		print "-- Fetching FASTA for ", $PFAM_ID, "\n";
		my @fasta = getFASTA($PFAM_ID);
		my @gene_names = humanGenesNames(@fasta);
		print "-- Finding GO terms por human genes\n";
		my %humanGOterms = GOterms($anno_mouse_file, $obo_file, $gene_names[0]);
		print "-- Finding similar genes in mouse annotation file\n";
		my %mouseGOgenes = findSimilarGenes($anno_mouse_file, %humanGOterms);
        my %frequent = moreSimilarGenes(%mouseGOgenes);
        
        my $filename = "mouseSim_results.txt";
        my $num = 30;
        print "-- First $num more similar mouse genes will be selected\n";
        print "-- Writing results in file $filename \n";
        
        open(my $fh, '>', $filename) or die "Could not open file '$filename' $!";
			print $fh "##############################################\n";
			print $fh "##                                          ##\n";
			print $fh "##              M o u s e S i m             ##\n";
			print $fh "##                                          ##\n";
			print $fh "##############################################\n";
			print $fh "\n\n";
            print $fh "\t GO TERMS \t GENE NAME\n";
            print $fh "\t----------\t----------\n";
            
        # Printing more frequent genes
        my $aux = 0;
        foreach my $gene (sort { $frequent{$b} <=> $frequent{$a} } keys %frequent) {    # Descending order for hash values
            $aux++;
            print $fh "\t$frequent{$gene} \t\t $gene\n";
            if ($aux == $num) {
                last;
            }          
        }
		close $fh;
    }
    else {
        print "-- ERROR: PFAM_ID is not valid. Aborting program.\n"
    }

}