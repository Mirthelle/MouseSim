#!/usr/bin/perl -w 

use strict; 
use LWP::Simple;

use IO::File;

use GO::TermFinder;
use GO::AnnotationProvider::AnnotationParser;
use GO::OntologyProvider::OboParser;

use GO::Utils::File    qw (GenesFromFile);
use GO::Utils::General qw (CategorizeGenes);


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
    print "-- Fetching FASTA for ", $acc, "\n";
    
    # Writing file
    my $filename = $acc.".fasta";
    open(my $fh, '>', $filename) or die "Could not open file '$filename' $!";
    print $fh $entry;
    close $fh;
    return $filename;
    print "-- Writing FASTA file for ", $acc, "\n";
}

######################################################################
sub humanGenesNames {                                               ##
######################################################################
## Using a fasta file as input writes a text file with human        ##
## gene names.                                                      ## 
######################################################################
    my($fasta_file) = @_;
    
    print "-- Selecting human gene names from file ", $fasta_file, "\n";
    open(FASTA, $fasta_file);
    
    # Selecting headers correspondig to human genes and storing them in @headers.
    my @headers = ();
    foreach my $line (<FASTA>) {
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
    return $filename; 
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


##############################################
### PRUEBAS -- MAIN PROGRAM                 ##
##############################################

#my $acc = "PF00870";
#my $fasta = "PF00870.fasta";
#getFASTA($acc);

#my $gene_names = humanGenesNames($fasta);
#print "@gene_names";

my $gene_names = "names.txt";
my $obo_file = "go-basic.obo";
my $anno_human_file = "gene_association.goa_human";
my $anno_mouse_file = "gene_association_mouse.mgi";
#my @hola = GOterms($gene_names, $anno_file, $obo_file);
#print scalar @hola, "\n";
#print "@hola";
#print scalar @{$hola[0]}, " ", scalar @{$hola[1]}, " ", scalar @{$hola[2]};

my %GOterms = GOterms($anno_mouse_file, $obo_file, $gene_names);
#print keys %GOterms, "\n";
#print %GOterms;
print "\n\n\nGO TERMS IN HUMAN ####################";
foreach my $i (keys %GOterms){
    print "ASPECT $i, ANNOTATED GENES: \n", join("\t", @{$GOterms{$i}});
}


my %terms = findSimilarGenes($anno_mouse_file, %GOterms);

print "\n\n\n SIMILAR GENES IN MOUSE ##################";
foreach my $i (keys %terms){
    print "############## ASPECT $i ################";
    foreach my $j (keys %{$terms{$i}}){
        print "$j =", join("\t", @{$terms{$i}{$j}}), "\n";
    }
    
}


#foreach my $i (keys %terms){
#    foreach my $j (keys %{$terms{$i}}){
#        print "\n", $j, " = ", join(" ", @{$terms{$i}{$j}});
#    }
#}

#foreach my $key (keys %{$terms{'C'}}){
#    print "\n", $key, "==", ${$terms{'C'}};
#}