=head1 NAME

SortCoordinates.pm - sort and merge arrays of regions

=head1 VERSION

version 0.1

=head1 SYNOPSIS

 use SortCoordinates;
 
 open (my $BED, "regions.bed");
 my @bed = <$BED>;

 my @sorted = SortCoordinates::sortByCoordinate
 (
     array => \@bed, 
     format => 'bed',
 );

 my @merged = SortCoordinates::mergeByCoordinate
 (
     array => \@sorted, 
     format => 'bed',
 );

=cut

package SortCoordinates;
use strict;
use warnings;
use Carp;

=head1 FUNCTIONS

=over 8

=item B<sortByCoordinate>

Given a reference to an array of regions this function sorts the array by contig and position.

Arguments

=over 12

=item array

A reference to an array of regions. 
=item format

The format of the regions. Default = BED. Valid formats are BED, VCF or regions (e.g. chr1:1000-2000).

=item value

The value for the given ID. If the INFO field is a FLAG this value should be ommited. 

=item contigs

Optional ARRAY or HASH reference of contig orders. Default behaviour is to sort contigs in ascibetical order.

If an ARRAY reference is passed, contigs will be sorted in the same order as this array. If a HASH reference is passed, it is expected that keys are contig names and values are integers representing the order of the contig.

=back

 my @sorted = SortCoordinates::sortByCoordinate(array => \@bed).
 
 my @sorted = SortCoordinates::sortByCoordinate
 (
     array => \@regions, 
     format => 'region'
 )'
 
 my @contigs = (1..22), qw /X Y M/;
 my @sorted = SortCoordinates::sortByCoordinate
 (
     array => \@bed, 
     contigs => \@contigs
 );

=cut

sub sortByCoordinate{
    my %args = @_;
    if (not $args{array}){
        croak "array argument is required for sortByCoordinate method ";
    }elsif( ref $args{array} ne 'ARRAY'){
        croak "array argument passed to sortByCoordinate must be an ARRAY ".
        "reference ";
    }
    
    my %contig_order = ();
    if ($args{contigs}){#optional order for contigs
        if (ref $args{contigs} eq 'ARRAY'){
            my $n = 0;
            %contig_order = map { $_ => $n++ } @{$args{contigs}};
        }elsif (ref $args{contigs} eq 'HASH'){
            %contig_order = %{$args{contigs}};
        }else{
            croak "contigs argument passed to sortByCoordinate must be either ".
           " an ARRAY or a HASH reference ";
        }
    }

    if (not $args{format} or uc($args{format}) eq 'BED'){
        return _sortBed( $args{array}, \%contig_order) ;
    }elsif(uc($args{format}) eq 'VCF'){
        return _sortVcf( $args{array}, \%contig_order) ;
    }elsif(uc($args{format}) eq 'REGION'){
        return _sortRegion( $args{array}, \%contig_order) ;
    }else{
        croak "unrecognised format ('$args{format}') ";
    }
}

sub _sortBed{
    my $array = shift;
    my $contigs = shift;
    return _sortArray
    (
        array     => $array,
        delimiter => "\t",
        seq       => 0,
        start     => 1,
        end       => 2,
        contigs   => $contigs,
    );
}

sub _sortVcf{
    my $array = shift;
    my $contigs = shift;
    return _sortArray
    (
        array     => $array,
        delimiter => "\t",
        seq       => 0,
        start     => 1,
        contigs   => $contigs,
    );
}

sub _sortRegion{
    my $array = shift;
    my $contigs = shift;
    return _sortArray
    (
        array     => $array,
        delimiter => '[:\-\t]',
        seq       => 0,
        start     => 1,
        end       => 2,
        contigs   => $contigs,
    );
}

sub _sortArray{
    my %args = @_;
    if ( %{$args{contigs}}){
        return sort { _byCoordinateWithContigs($a, $b, %args) } @{$args{array}};
    }else{
        return sort { _byCoordinate($a, $b, %args) } @{$args{array}};
    }
}

sub _byCoordinateWithContigs{
    my ($a, $b, %args) = @_;
    chomp $a;
    chomp $b;
    my @a_split = split(/$args{delimiter}/, $a);
    my @b_split = split(/$args{delimiter}/, $b);
    my $a_chrom = $args{contigs}->{$a_split[$args{contig}]}; 
    my $b_chrom = $args{contigs}->{$b_split[$args{contig}]};
    
    if (not defined $a_chrom or not defined $b_chrom){
    #if either chrom is not in our supplied contigs put them 
    #to the end of the array
        if (defined $a_chrom){
            return -1;
        }elsif(defined $b_chrom){
            return 1;
        }else{
        #if both chrom are not in contigs give them arbitrary number 
        #and continue to sort on coordinate
            $a_chrom = 0;
            $b_chrom = 0;
        }
    } 

    if ($args{end}){
        return 
        (
            $a_chrom               <=> $b_chrom               ||
            $a_split[$args{start}] <=> $b_split[$args{start}] ||
            $a_split[$args{end}]   <=> $b_split[$args{end}]   ||
            $a                     cmp $b
        );
    }else{
        return 
        (
            $a_chrom               <=> $b_chrom               ||
            $a_split[$args{start}] <=> $b_split[$args{start}] ||
            $a                     cmp $b
        );
    }
}

sub _byCoordinate{
    my ($a, $b, %args) = @_;
    my @a_split = split(/$args{delimiter}/, $a);
    my @b_split = split(/$args{delimiter}/, $b);
    if ($args{end}){
        return 
        (
            $a_split[$args{seq}]   cmp $b_split[$args{seq}]   ||
            $a_split[$args{start}] <=> $b_split[$args{start}] ||
            $a_split[$args{end}]   <=> $b_split[$args{end}]   ||
            $a                     cmp $b
        );
    }else{
        return 
        (
            $a_split[$args{seq}]   cmp $b_split[$args{seq}]   ||
            $a_split[$args{start}] <=> $b_split[$args{start}] ||
            $a                     cmp $b
        );
    }
}
=item B<mergeByCoordinate>

Given a reference to a sorted array of regions this function merges the array by contig and position.

Arguments

=over 12

=item array

A reference to an array of regions. This must already be sorted in order for the merge to work properly.

=item format

The format of the regions. Default = BED. Valid formats are BED or regions (e.g. chr1:1000-2000).

=back

 my @merged = mergeByCoordinate(array => \@sorted).

=cut 


sub mergeByCoordinate{
    my %args = @_;
    if (not $args{array}){
        croak "array argument is required!\n";
    }
    if (not $args{format} or uc($args{format}) eq 'BED'){
        return _mergeBed( $args{array}, $args{keep_info} ) ;
    }elsif(uc($args{format}) eq 'REGION'){
        return _mergeRegion( $args{array}, $args{keep_info} ) ;
    }else{
        croak "unrecognised format ('$args{format}') ";
    }
}

sub _mergeBed{
    my $array = shift;
    my $keep_info = shift;
    return _mergeArray
    (
        array     => $array,
        delimiter => "\t",
        seq       => 0,
        start     => 1,
        end       => 2,
        keep_info => $keep_info,
    );
}
sub _mergeRegion{
    my $array = shift;
    my $keep_info = shift;
    return _mergeArray
    (
        array     => $array,
        delimiter => '[:\-\t]',
        seq       => 0,
        start     => 1,
        end       => 2,
        keep_info => $keep_info,
    );
}

sub _mergeArray{
    my %args = @_;
    my %prev_region = ();
    my @merged = ();
    foreach my $reg (@{$args{array}}){
        chomp $reg;
        my @split = split(/$args{delimiter}/, $reg);
        if (not %prev_region){
            %prev_region = _convertRegionToHash(\@split, %args);
        }else{
            if 
            (
                $split[$args{seq}] eq $prev_region{seq} and
                $split[$args{start}] <= $prev_region{end} and
                $split[$args{end}] >= $prev_region{start} 
                #the final comparison shouldn't be necessary if array is sorted
                #but is included so as not to completely FUBAR any unsorted 
                #array passed to this method
            ){#merge overlapping regions
                if ($split[$args{end}] > $prev_region{end}){
                    $prev_region{end} = $split[$args{end}];
                }
                $prev_region{extra} .= "," . join("|", @split[$args{end}+1..$#split]);
            }else{#not overlapping
                push @merged, _getRegion(\%prev_region, %args);
                %prev_region = _convertRegionToHash(\@split, %args);
            }
        }
    }
    push @merged, _getRegion(\%prev_region, %args);
    return @merged;
}

sub _getRegion{
    my ($reg, %args) = @_;
    my $line = '';
    if ($args{delimiter} eq '[:\-\t]'){
        $line = "$reg->{seq}:$reg->{start}-$reg->{end}";
    }else{
        $line = join("\t", map { $reg->{$_} } qw/ seq start end /);
    }
    if ($args{keep_info}){
        $line .= "\t$reg->{extra}";
    }
    return $line;
}

sub _convertRegionToHash{
    my ($split, %args) = @_;
    my %region = ();
    foreach my $k (qw /seq start end /){
        $region{$k} = $split->[$args{$k}];
    }
    $region{extra} = join("|", @$split[$args{end}+1..$#{$split}]);
    return %region;
}

=back

=head1 AUTHOR

David A. Parry

=head1 COPYRIGHT AND LICENSE

Copyright 2016  David A. Parry

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version. This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details. You should have received a copy of the GNU General Public License along with this program. If not, see <http://www.gnu.org/licenses/>.

=cut

1;


