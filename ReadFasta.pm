package ReadFasta;
use strict;
use warnings;
use Exporter;
our $VERSION = 0.1;
#use vars qw($VERSION @ISA @EXPORT_OK );

our @ISA = qw(Exporter);
our @EXPORT_OK = qw( read_fasta );

##################################################
sub read_fasta{
#returns a hash of sequence IDs to sequence
    my $f = shift;
    open (my $FH, $f) or die "Could not read $f: $!\n";
    my %fasta = (); 
    my %buffer = ();
    while (my $l = <$FH>){
        chomp $l; 
        next if not $l;
        if ($l =~ /^>/){
            if (%buffer){
                $fasta{$buffer{name}} = $buffer{seq};
                %buffer = (seq => '');
            }
            $l =~ s/>//;
            $buffer{name} = $l;
        }else{
            $buffer{seq} .= $l;
        }
    }
    if (%buffer){
        $fasta{$buffer{name}} = $buffer{seq};
    }
    close $FH;
    return %fasta;
}

1;

=head1 NAME

ReadFasta.pm - read a fasta file and return a hash of sequence IDs to names

=head1 SYNOPSIS

    use ReadFasta qw ( read_fasta ); 

    my %fastqs = read_fasta("input.fasta");
    foreach my $k (keys %fastqs){
        print "For ID $k sequence is $fastqs{$k}\n";
    }

=head1 DESCRIPTION

A simple, lightweight module to read in a whole fasta file and return a hash of 
sequence IDs to sequences.

Designed for convenience when handling relatively small sets of sequences. Not 
recommended for large files.

=head1 AUTHOR

David A. Parry

=head1 COPYRIGHT AND LICENSE

Copyright 2017 David A. Parry

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version. This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details. You should have received a copy of the GNU General Public License along with this program. If not, see <http://www.gnu.org/licenses/>.

=cut


