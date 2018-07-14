#! /usr/bin/env perl
use Cwd;
use Getopt::Long;
use File::Basename;
use strict;
use warnings;

$ENV{PATH} = "/usr/local/bin:/bin:/usr/bin:/usr/local/sbin:/usr/sbin:/sbin:/usr/kerberos/bin";

# Some defines for our returns
my $ALLSWELL = 0;
my $FAILEDCURLWRAPPER = 1;
my $FAILEDCURLRUNTIME = 2;
my $CHECKSUMFAILED = 3;
my $NORUNTIME = 4;
my $RINSTALLFAIL = 5;
my $UNTARRUNTIME = 6;
my $UNTARJOBDATA = 7;
my $FAILEDRLIBUNTAR = 8;
my $FAILEDRUN = 9;
my $MISSINGOUTPUT = 10;
my $BADUSAGE = 11;
my $NOCHECKSUM = 12;
my $FAILEDCURLCHECKSUM = 13;
my $FAILEDCURLSLIBS = 14;

my $unique = shift @ARGV;
my $location = getcwd();
my $Rlocation = "$location/R";
my $renviron;
my $dbinstalllog =  "runR.$unique.out";
print "Trying to open logfile... $dbinstalllog\n";
#open(OLDOUT, ">&STDOUT");
#open(OLDERR, ">&STDERR");
open(STDOUT, ">>$dbinstalllog") or die "Could not open $dbinstalllog: $!";
open(STDERR, ">&STDOUT");
select(STDERR);
 $| = 1;
 select(STDOUT);
  $| = 1;


my $stampfilename = "AuditLog.$unique";
my $cmdtorun = shift @ARGV;
my $currentslibs = "SLIBS1.tar.gz";

if(!(-f "$currentslibs")) {
	my $localcurlres = FetchViaCurl("$currentslibs","$unique.out","same");
	if($localcurlres != 0) {
		exit($FAILEDCURLSLIBS);
	}
	system("tar -zxvf $currentslibs"); 
}




my $prebuilt = 'built-sl5-R-2.10.1.tar.gz';
my $installfrom = 'R-2.10.1';
my $curlres = FetchViaCurl($prebuilt, "$unique.out", "same");

        if($curlres != 0) {
                    system("rm $prebuilt");
                            print "Fetch of runtime <$prebuilt> failed via curl\n";
                                    exit($FAILEDCURLRUNTIME);
                                        }
system("tar -zxvf ./$prebuilt");
system("chmod +rwx $prebuilt");
#system("chmod -R a+x .");



# do we have prebuilt libraries?
if(-f "SLIBS.tar.gz") {
	my $libres = system("tar -zxvf SLIBS.tar.gz");
	if($libres != 0) {
		return($FAILEDRLIBUNTAR);
	}
}
# This might have come from SLIBS.tar.gz (immediately above)
# or $currentslibs("SLIBS1.tar.gz") around line 134.
if( -d "$location/SS" ) {
	$ENV{LD_LIBRARY_PATH} = "$location/SS:$ENV{LD_LIBRARY_PATH}";
}

# Check for gcc
my $gcc = system("which gcc");
if ($gcc != 0) {
	print "no gcc available; R install failure\n";
	chdir("$location");
	system("rm -r R-2.10.1 SS AuditLog.* built-sl5-R-2.10.1.tar.gz SLIBS1.tar.gz");
	exit($RINSTALLFAIL);
}
chdir("$location");

# R needs to prepare runtime
chdir("$installfrom");
my $Rres = system("make prefix=$Rlocation install");

chdir("$location");

# Set up R environment
$ENV{R_HOME} = "$Rlocation";
$ENV{LD_LIBRARY_PATH} = "$Rlocation/lib:$ENV{LD_LIBRARY_PATH}";

# do we have prebuilt libraries?
if(-f "RLIBS.tar.gz") {
	my $libres = system("tar -zxvf RLIBS.tar.gz");
	if($libres != 0) {
		return($FAILEDRLIBUNTAR);
	}
	$ENV{HOME} = $location;
	$renviron = "$location/.Renviron";
	open(RLIB,">$renviron") or die "Can not create <$renviron>:$!\n";
	print RLIB "R_LIBS_USER=$location/RR/library\n";
	close(RLIB);
}



my $runthis = "$Rlocation/bin/Rscript --no-save ./$cmdtorun @ARGV";
print "about to execute <<$runthis>>\n";
my $res = system("$runthis");
if($res != 0) {
	print "Execution of <$runthis> failed\n";
} else {
	print "Done Running. Result was <$res>\n";	

}

system("rm AuditLog.* built-sl5-R-2.10.1.tar.gz SLIBS1.tar.gz");

exit(0);

sub write_audit_message {
        open(STAMP,">>$stampfilename");
            my($type, %classad) = @_;
                print STAMP ( timestamp(), " $type ");
                    foreach my $key (sort keys (%classad)) {
                                my $filter = $classad{$key};
                                        $filter =~ s/>/GT/g;
                                                print STAMP "$key=<$filter> ";
                                                    }
                                                        print STAMP "\n";
                                                            close(STAMP);
}

sub timestamp {
        return scalar localtime();
}


sub FetchViaCurl
{
	my $target = shift;
	my $unique = shift;
	my $savename = shift;

	my $tries = 30;
	my $trycount = 0;
	my $sleeptime = 6;
	my $starttime = time();
	my $random = time();
	my $url = "";
		$url = "http://proxy.chtc.wisc.edu/SQUID/";
	print "URL: <<$url>>\n";

	while($trycount < $tries) {
		my $curlcmd = "";
		if($savename eq "same") {
			$curlcmd = "curl -H \"Pragma:\" -o $target $url$target";
		} else {
			$curlcmd = "curl -H \"Pragma:\" -o $savename $url$target";
		}
		print "Curl attempt:<$curlcmd>\n";
	    $curlres = system("$curlcmd");
		# trace data per attempt
		my $stoptime = time();
		my $elapsed = $stoptime - $starttime;
		my $fetchfile = "$url$target";
		my $result = $curlres >>=8;
		my $signal = $curlres & 255;
		write_audit_message("CURL",(url=>$fetchfile,time=>$elapsed,signal=>$signal,result=>$result));
		if($curlres != 0) {
			$trycount += 1;
			sleep($sleeptime);
		} else {
			last;
		}
	}
	print "Curl error <$curlres>\n";
	return($curlres);
}
