#!/usr/bin/perl -w
 
$| = 1;

# Program name
print "\n--- $0 ---\n";

# Collect command line arguments
$USAGE = "perl $0 <clusters.fasta> <gap_size(default:10000)>\n";
unless ($ARGV[0]) {
	print "\nUsage: $USAGE\n";
	exit;
}
$infile = $ARGV[0];
if ($ARGV[1]) {
	$max_gap = $ARGV[1];
} else {
	$max_gap = 10000;
}
$outfile = $infile.'.merge.fas';

# Processing message
print "Processing";

open (IN, $infile) or die "\nCannot open file \'$infile\'!\n";
while ($line = <IN>) {
	$line =~ s/\s+$//; #chomp
	
	if ($line =~ /directionality/) {
	
		$cluster++;
		
		@d = split('\s+', $line);
		
		$chr = $d[0];
		$loc = $d[1];
		$dir = $d[5];
		
		@loc_crd = split('-', $loc);
		$loc_beg = $loc_crd[0];
		$loc_end = $loc_crd[1];
		
		$cluster_chr{$cluster} = $chr;
		$cluster_loc_beg{$cluster} = $loc_beg;
		$cluster_loc_end{$cluster} = $loc_end;
		$cluster_dir{$cluster} = $dir;
		
	} else {
		
		$cluster_seq{$cluster} = $line;
		
	}
	
}
close(IN);

# Search for merge points
foreach my $clu (sort { $a <=> $b } keys %cluster_seq) {
	
	$clu_next = $clu+1;
	
	if ($cluster_chr{$clu_next}) {
		if ($cluster_chr{$clu} eq $cluster_chr{$clu_next}) {
		
			$gap = $cluster_loc_beg{$clu_next}-$cluster_loc_end{$clu};
		
			if ($gap <= $max_gap) {
				
				$merge_point = ($clu+0.5);
				$mps{$merge_point} = 1;
				$mp_gaps{$merge_point} = $gap;
				
			}
		}
	}
}

open (OUT, ">$outfile") or die "\nCannot open file \'$outfile\'!\n";

# Go through every cluster
foreach my $clu (sort { $a <=> $b } keys %cluster_seq) {
	
	# Check if merge point exists at cluster
	$mp = ($clu+0.5);
	if ($mps{$mp}) {
	
		# skip already considered adjacent merge points
		unless ($adj_mps{$mp}) {
		
			# Add cluster to list for skipping
			$proc_clus{($mp+0.5)} = 1;
		
			# Merge two clusters at merge point
			$clu_merged = ($mp-0.5).'_'.($mp+0.5);
			$Nfill = 'N'x$mp_gaps{$mp};
			$clu_seq_merge = $cluster_seq{($mp-0.5)}.$Nfill.$cluster_seq{($mp+0.5)};
			
			# Go to putative adjacent merge point
			$adj_mp = $mp+1;
			
			# Check if adjacent merge point exists
			while ($mps{$adj_mp}) {
			
				# Add cluster to list for skipping
				$proc_clus{($adj_mp+0.5)} = 1;
			
				# if adjacent merge point exists 
				# merge next cluster with previously merged clusters
				$clu_merged = $clu_merged.'_'.($adj_mp+0.5);
				$Nfill = 'N'x$mp_gaps{$adj_mp};
				$clu_seq_merge = $clu_seq_merge.$Nfill.$cluster_seq{($adj_mp+0.5)};
				
				# Add adjacent merge point to list for skipping
				$adj_mps{$adj_mp} = 1;
				
				# Go to next putative adjacent merge point
				$adj_mp = $adj_mp+1;
			}
			
			# Save merged cluster and info
			$clu_m_fir = $clu_merged;
			$clu_m_fir =~ s/_.*//;
			$clu_m_las = $clu_merged;
			$clu_m_las =~ s/.*_//;
			print OUT ">$clu_merged\t$cluster_chr{$clu_m_fir}\t$cluster_loc_beg{$clu_m_fir}-$cluster_loc_end{$clu_m_las}\t$cluster_dir{$clu_m_fir}\/$cluster_dir{$clu_m_las}\n$clu_seq_merge\n";
		}
	} else {
	
		# Save non-merged clusters and skip merged clusters
		unless ($proc_clus{$clu}) {
			print OUT ">$clu\t$cluster_chr{$clu}\t$cluster_loc_beg{$clu}-$cluster_loc_end{$clu}\t$cluster_dir{$clu}\n$cluster_seq{$clu}\n";
		}
	}
}
close(OUT);

print " done.\n";

exit;