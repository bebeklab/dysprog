use warnings;
use strict;
use Array::Shuffle qw (shuffle_array);
use Getopt::Long;
use Data::Dumper;
$|++;


###################################################################################################################################
sub pr_err{
    print "Usage:\n",
    "\t--fsm_dir <FSM Output Directory>\n",
    "\t--expr_file <Log Scale Expression File>\n",
    "\t--clin_file <Clinical File>\n",
    "\t--output <Base Name for Output>\n",
    "\t--thr <Maximum Match Count for Digraph Comparison>\n";
}

sub mean{
    my $arrayRef = shift;
    my $total = 0;
    my $count = 0;
    foreach (@{$arrayRef}){
        if(defined $_ && $_ ne '' && $_ ne 'NA'){
            $total+=$_;
            $count++;
        }
    }
    return ($total/$count);
}

sub stdev{
    my $arrayRef = shift;
    my $mean = mean($arrayRef);
    my $score = 0;
    my $count = 0;
    foreach (@{$arrayRef}){
        if(defined $_ && $_ ne '' && $_ ne 'NA'){
            $score+=($_ - $mean)**2;
            $count++;
        }
    }
    return (sqrt($score/($count-1)));
}

sub zscore{
    my $arrayRef = shift;
    my ($mean, $sd) = (mean($arrayRef), stdev($arrayRef));
    my @zscores = map{(defined $_ && $_ ne '' && $_ ne 'NA')?(($_ - $mean)/$sd):('NA')}(@{$arrayRef});
    return \@zscores;
}
sub z_transform{
    my $expr_data = shift;
    my @patients = keys %{$expr_data};
    my %genes = ();
    foreach my $pat (@patients){
        map{$genes{$_} = 1}(keys %{$expr_data->{$pat}});
    }
    my %expr_standard = ();
    foreach my $gene (keys %genes){
        my @values = ();
        map{push @values, $expr_data->{$_}{$gene}}(@patients);
        my $zscores = zscore(\@values);
        map{$expr_standard{$patients[$_]}{$gene} = $zscores->[$_]}(0..@patients-1);   
    }
    return \%expr_standard;
}
sub parseExpr{
    my $file = shift;
    my (%Expr, @patients),
    open my $fH, "<$file" || die "Error! Failed to read file: $file\t$!\n";
    while(<$fH>)
    {
        chomp;
        if($. == 1){
            @patients = split /\t/;
            next;
        }
        my @line = split /\t/;
        my $gene = shift @line;
        map{$Expr{$patients[$_]}{$gene} = $line[$_]}(0..@line-1);
    }
    close $fH;   
    return \%Expr;
}

sub choose_best{
    my ($graph_gene_data, $graph_data, $size_data) = @_;
    my @sorted_runs = sort{$size_data->{$b} <=> $size_data->{$a}}(keys %{$size_data});
    
    my $max = $size_data->{$sorted_runs[0]};
    my @max_runs = grep{$size_data->{$_} == $max}(keys %{$size_data});
    if(@max_runs == 1){
        return $max_runs[0];
    }else{
        #print "\nCalculating Matching Fractions\n";
        my %fract_data = ();
        foreach my $run (@max_runs){
            my @graphs = split /,/, $graph_data->{$run};
            my %processed = ();
            my $total = 0;
            foreach my $dig1 (@graphs){
                foreach my $dig2 (@graphs){
                    next if $dig1 eq $dig2;
                    next if exists $processed{"$dig1-$dig2"} || exists $processed{"$dig2-$dig1"};
                    my $match = grep{exists $graph_gene_data->{$dig2}{$_}}(keys %{$graph_gene_data->{$dig1}});
                    my $frac = $match / (scalar(keys %{$graph_gene_data->{$dig1}}) + scalar(keys %{$graph_gene_data->{$dig1}}));
                    $total+=$frac;
                    $processed{"$dig1-$dig2"} = 1;
                }
            }
            $fract_data{$run} = $total;
        }
        #map{print $_, "\t", $fract_data{$_}, "\n"}(keys %fract_data);
        my @sorted_fract_data = sort{$fract_data{$a} <=> $fract_data{$b}}(keys %fract_data); # Choose Smallest Match
        return $sorted_fract_data[0];
    }
}

sub scoreDigraph{
    print " % Calculating Digraph Scores\n";
    my ($Expr, $Digraph_Data) = @_;
    my %DigraphScores = ();
    my $total_patients = scalar(keys %{$Expr});
    my $patient_count = 0;
    foreach my $pat (keys %{$Expr}){
        foreach my $path (keys %{$Digraph_Data}){
            my @ints = keys %{$Digraph_Data->{$path}};
            my ($score, $count) = (0, 0);
            map{
                    my ($gene1, $gene2) = split />/;
                    if(exists $Expr->{$pat}{$gene1} && exists $Expr->{$pat}{$gene2} &&
                    $Expr->{$pat}{$gene1} ne "NA" && $Expr->{$pat}{$gene2} ne "NA"){
                    $score += sqrt($Expr->{$pat}{$gene1}**2 + $Expr->{$pat}{$gene2}**2);
                    $count++;
                    }
                    
                }(@ints);
            ($count > 0)?($DigraphScores{$path}{$pat} = $score / $count):($DigraphScores{$path}{$pat} = 'NA');
        }
        $patient_count++;
        my $percent = $patient_count*100/$total_patients;
        printf "\r\e[K   --> Status: %.2f", $percent if($patient_count % 10 == 0);
    }
    return \%DigraphScores;
}
###################################################################################################################################


my ($fsm_rep_dir, $expr_file, $outfile, $clin_file, $thr);

GetOptions("fsm_dir=s" => \$fsm_rep_dir,
           "expr_file=s" => \$expr_file,
           "output=s" => \$outfile,
           "clin_file=s" => \$clin_file,
           "thr=i" => \$thr);


if(!defined $fsm_rep_dir || $fsm_rep_dir eq '' || !-e $fsm_rep_dir || !-d $fsm_rep_dir){
    print "!Error in FSM Result Directory\n";
    pr_err;
    exit 1;
}elsif(!defined $expr_file || $expr_file eq '' || !-e $expr_file){
    print "!Error in Expression File\n";
    pr_err;
    exit 1;
}elsif(!defined $thr || $thr eq '' || $thr >= 8){
    print "!Error in threshold level (Not defined or > 8)\n";
    pr_err;
    exit 1;
}elsif(!defined $clin_file || $clin_file eq '' || !-e $clin_file){
    print "!Error in Clinical File\n";
    pr_err;
    exit 1;
}

printf " %% Reading FSM Output Directory\n";
opendir(my $dirH, $fsm_rep_dir) || die $!;
my @local_dirs = grep{m/max8.+?\d+\.fsgs$/}(readdir $dirH);
closedir $dirH;
if(@local_dirs == 0){
    print "!No sub-directories found\n";
    exit 1;
}
shuffle_array @local_dirs;
my $nameCount = 0;
my %digraph_data = ();
my %digraph_gene_data = ();
my %genes = ();
foreach my $dir_idx (0..@local_dirs-1){
    my $dir = $local_dirs[$dir_idx];
    print "   --> Reading $dir\n";
    opendir $dirH, $fsm_rep_dir.'/'.$dir || die $!;
    my @files = grep{m/patterns\.dot/}(readdir $dirH);
    closedir $dirH;
    if(@files != 1){
        print "!Directory $fsm_rep_dir/$dir should contain a single patterns.dot file, skipping\n";
        next;
    }
    open my $fH, "<", $fsm_rep_dir.'/'.$dir.'/'.$files[0] or die $!;
    my $file = "";
    $file.=$_ while <$fH>;
    close $fH;
    while($file=~m/digraph.+?\{(.+?)\}/sg){
        my $graph="Digraph-$nameCount";
        my $localDigraph = $1;
        my %labels = ();
        while($localDigraph=~m/\[idx=\"(.+?)\",label=\"(.+?)\"\]/sg){
            $labels{$1} = $2;
            $genes{$2} = 1;
        }
        while($localDigraph =~ m/(\d+)\s+->\s+(\d+)/sg){
            my ($prot1, $prot2) = ($labels{$1}, $labels{$2});
            next if !defined $prot1 || !defined $prot2;
            $digraph_data{$graph}{$prot1.">".$prot2} = 1;
            $digraph_gene_data{$graph}{$prot1} = 1;
            $digraph_gene_data{$graph}{$prot2} = 1;
        }
        $nameCount++;
    }
}

printf " %% Parsed %d Digraphs, %d Genes\n", scalar(keys %digraph_data), scalar(keys %genes);

print " % Searching for maximum set\n ";
my %all_results = ();
my %all_results_size = ();
my $total_gene_thr = scalar(keys %genes) * 0.95;
RUN: foreach my $run (0..1000){
	my @graphList = keys %digraph_data;
	my %filtered_graph_list = ();
	my %local_genes = ();
	while(scalar(keys %local_genes) < $total_gene_thr && scalar(@graphList) > 0){
		my $idx = int(rand(scalar(@graphList)));
		my $dig = splice @graphList, $idx, 1;
		my $match = grep{exists $local_genes{$_}}(keys %{$digraph_gene_data{$dig}});
		next if $match > $thr;
		$filtered_graph_list{$dig} = 1;
		map{$local_genes{$_} = 1}(keys %{$digraph_gene_data{$dig}});
                print "\r\e[K   --> Run:$run, Local Size:", scalar(keys %filtered_graph_list), ", Local Gene Size:", scalar(keys %local_genes);
	}
    $all_results{sprintf("Run%d", $run)} = join(",", keys %filtered_graph_list);
    $all_results_size{sprintf("Run%d", $run)} = scalar(keys %filtered_graph_list);
    print "\n";
}

my $best_run = choose_best(\%digraph_gene_data,\%all_results,\%all_results_size);
print " [ Best Run: $best_run - ", $all_results_size{$best_run}," Digraphs ]\n";

## Extract Digraph Data for Best Set ##
my %Digraph_Data = ();
my @digraph_list = split /,/, $all_results{$best_run};
foreach my $dig (@digraph_list){
    map{$Digraph_Data{$dig}{$_} = 1}(keys %{$digraph_data{$dig}});
}

print " % Reading Expression File ";
my $Expr = parseExpr($expr_file);
my $expr_data = z_transform($Expr);
printf "[ %d Patients, %d Genes ]\n", scalar(keys(%{$Expr})), scalar(keys(%{$Expr->{(keys(%{$Expr}))[0]}}));


my $Digraph_Scores = scoreDigraph($expr_data, \%Digraph_Data);
my $Dcount = scalar(keys(%{$Digraph_Scores}));
my ($scoreF, $digF) = ($outfile."-Scores_D$Dcount", $outfile."-Digraphs_D$Dcount");
print "\n";



print " % Writing Results";
my @patients = keys %{$expr_data};
open my $fH, ">", $scoreF || die "!Cannot Write File $scoreF\n";
print $fH join("\t", @patients), "\n";
foreach my $path (keys %{$Digraph_Scores}){
    print $fH $path;
    map{(exists $Digraph_Scores->{$path}{$_} && $Digraph_Scores->{$path}{$_} ne "NA")?(printf $fH "\t%.3f", $Digraph_Scores->{$path}{$_}):(print $fH "\tNA")}(@patients);
    print $fH "\n";
}
close $fH;

open $fH, ">", $digF || die "!Cannot Write File $digF\n";
foreach my $graph (keys %{$Digraph_Scores}){
    print $fH $graph, "\t", join(',', keys(%{$Digraph_Data{$graph}})), "\n";
}
close $fH;
print "\n";

exit 0; ## Temporary

if(!-e "nmf_cluster.R"){
    print "!Could not locate \'nmf_cluster.R\' script\n";
    exit 1;
}

print " % Running NMF Clustering\n";
my @cluster_args = ("Rscript", "~/FSM-Run/nmf_cluster.R", $scoreF, $clin_file, $outfile, "1500");
exec @cluster_args;

exit 0;





