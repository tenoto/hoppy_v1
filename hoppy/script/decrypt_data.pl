#!/usr/bin/perl

#############################################################################
#   Description:                                                            #
#                                                                           #
#   This script decrypts all *.pgp or *.gpg files in a directory tree       #
#   or individual files.                                                    #
#                                                                           #  
#   The script tests if you have PGP and/or GPG installed.                  #
#   If the files end in .gpg, GPG will be used to decrypt.                  #
#   If only GPG is available, GPG will be used to decrypt.                  #
#   If the -gpg option is used, GPG will be used to decrypt.                #
#   Otherwise PGP will be used to decrypt.                                  #
#                                                                           #
#   The script has to work at least in these platforms and with these       #
#   versions of PGP and GPG:                                                #
#                                                                           #
#   For Linux and Sun, it works with PGP Version 6.5.8.(Command="pgp")      #
#                                                                           #
#   For Linux, it works with GPG (GnuPG) Version 1.0.7 (Command="gpg")      # 
#                                                                           #
#   For Macintosh, it works with GPG (GnuPG) Version 1.2.4 (Command="gpg")  #
#                                    (GnuPG/MacGPG2) 2.0.22                 #
#                  (Supported algorithms:                                   #
#                  Pubkey: RSA, RSA-E, RSA-S, ELG-E, DSA, ELG               #
#                  Cipher: IDEA, 3DES, CAST5, BLOWFISH, AES, AES192,        #
#                          AES256, TWOFISH                                  #
#                  Hash:   MD5, SHA1, RIPEMD160, SHA256).                   #
#                                                                           #
#                                                                           #
#   The extension of the files to be decrypted has to be .pgp or .gpg       #
#   2014/Jul Update version to allow Mac GPG                                #
#############################################################################

use strict;
use Getopt::Long;
use File::Find;

#my $VERSION = "3.2.1";
my $VERSION = "3.3.1";

my $print_version = 0;
my $print_help = 0;
my $remove = 0;
my $password = undef;
my $directory = undef;
my $gpg = 0;
my $forcepgp = 0;
my $quiet = 0;
GetOptions ("version"     => \$print_version,
            "help"        => \$print_help,
            "remove"      => \$remove,
            "password=s"  => \$password,
            "directory=s" => \$directory,
            "forcepgp"    => \$forcepgp,
            "gpg"         => \$gpg,
            "quiet"       => \$quiet );

######################################
#           Print Help               #
######################################

if($print_help){

       &PRINT_HELP;
       exit;
}

######################################
#           Print Version            #
######################################

if($print_version){

       print"decrypt.pl version $VERSION\n";
       exit;
}
    
######################################
#     Determine files to decrypt     #
######################################

if ( !defined $directory && @ARGV == 0 ) {
    print "Enter the Directory ('.' for current directory): \n";
    $directory=<STDIN>;
    chomp($directory);
}
if ( defined $directory && @ARGV ) {
   die "Must enter directory (-d option) or list of files, not both\n";
}
my @decrypt_files = ();
if ( @ARGV ) {
   foreach my $file (@ARGV) {
      if ( $file =~ /\.gpg$/ ) {
         $gpg = 1;
         push @decrypt_files, $file;
      } elsif ( $file =~ /\.pgp$/ ) {
         push @decrypt_files, $file;
      }
   }
} else {
   find({wanted => \&visit, follow => 1}, $directory);
}
if ( @decrypt_files == 0 ) {
   die "No files to decrypt\n";
}
######################################
#     Test PGP and GPG versions      #
######################################

my $pgpworks=0;
my $pgpversion = undef;

open PGPTEST, "pgp -h 2>&1 |";
while ( <PGPTEST> ) {
   chomp;
   if ( $_ =~ /Version (.*)$/ ) {
      $pgpworks = 1;
      $pgpversion = $1;
      last;
   }
}
close PGPTEST;

my $gpgworks=0;
my $gpgversion = undef;

open GPGTEST, "gpg --version 2>&1 |";
while ( <GPGTEST> ) {
   print "got gpgtest content  $_\n";
   chomp;
#
# new statement to support MAC version 2 
   if ( $_ =~ /^gpg.+\(GnuPG.*\) (.*)$/ ) {
#   if ( $_ =~ /^gpg \(GnuPG\) (.*)$/ ) {
      $gpgworks = 1;
      $gpgversion = $1;
   }
}
#print "force  $gpg \n";
close GPGTEST;

######################################
#     Determine decryption command   #
######################################
if ( !$pgpworks ) { $gpg = 1; }
#print " PGP TEST:   $gpgworks , ver $gpgversion , for $gpg\n" ;
if ( $forcepgp )  { $gpg = 0; }
#print " FORCE TEST:   $gpgworks , ver $gpgversion , for $gpg\n" ;
if ( $gpg && !$gpgworks )  { die "GPG is unavailable\n"; }
if ( $forcepgp && !$pgpworks )  { die "PGP is unavailable\n"; }

my $decrypt_cmd = undef;
if ( $gpg ) {
   $decrypt_cmd = "gpg --batch --no-tty --passphrase-fd 0";
   if ( !$quiet ) { print "\nUsing GPG $gpgversion\n\n"; }
} else {
   $ENV{'PGPPASSFD'} = '0';
   $decrypt_cmd = "pgp +batchmode";
   if ( !$quiet ) { print "\nUsing PGP $pgpversion\n\n"; }
}

######################################
#           DECRYPT DATA             #
######################################

	
foreach my $file (@decrypt_files) {

   my $newfile = undef;
   if ( $file =~ /^(.*)\.(gpg|pgp)$/ ) {
      $newfile = $1;
   } else {
      if ( !$quiet ) { print "WARNING: $file not encrypted, skipped...\n"; }
      next;
   }
   if ( -e $newfile ) {
      if ( !$quiet ) { 
         print "WARNING: $newfile already exists, skipped decrypt...\n";
      }
      next;
   }

   open DECRYPT, "| $decrypt_cmd -o $newfile $file 1>/dev/null 2>/dev/null" or 
                                        die "Failed to run $decrypt_cmd";
   if ( !defined $password ) {
      print "Enter the password: \n";
#      system('/bin/stty', '-echo');
      $password=<STDIN>;
      chomp($password);
#      system('/bin/stty', 'echo');
   }
   print DECRYPT "$password\n";
   close DECRYPT;
   my $retcode = $?;
   if ( $retcode == 0 || $retcode > 0x80 ) {
      if ( $retcode > 0x80 ) { $retcode >>= 8; }
      if ( $retcode != 0 && $retcode != 1 ) {
         die "Failed to decrypt $file\n";
      }
      if ( -e $newfile ) {
         if ( !$quiet ) { print "Decrypted: $newfile\n"; }
         if ( $remove ) { unlink $file; }
      }
   } elsif ( $retcode == 0xff00 ) {
      die "FAILED: $decrypt_cmd $file, error is $!\n";
   } else {
      die "FAILED: $decrypt_cmd $file, core dump: $retcode\n";
   }
}
######################################
#               HELP                 #
######################################

sub PRINT_HELP{
    
print "\n";
print "This script decrypts all *.pgp or *.gpg files in a directory tree\n";
print "specified with the -d option.  Individual files may be decrypted\n";
print "by giving them as arguments to the script.\n";
print "\n"; 
print "The script tests if you have PGP and GPG installed.\n";
print "By default the script uses PGP, however, for any of the following cases\n";
print "GPG is used:\n";
print "  * The -g option is set\n";
print "  * The files end in .gpg\n";
print "  * PGP is unavailable\n";
print "\n";                                                                          
print "The script has to work at least in these platforms and with these\n";
print "versions of PGP and GPG: \n";
print "\n";
print "For Linux and Sun, it works with PGP Version 6.5.8.(Command='pgp')\n";
print "For Linux, it works with GPG (GnuPG) Version 1.0.7 (Command='gpg')\n";
print "For Macintosh, it works with GPG (GnuPG) Version 1.2.4 (Command='gpg')\n";
print "                                  (GnuPG/MacGPG2) Vaersion 2.0.22  \n";
print "\n";
print "The extension of the files to be decrypted has to be .pgp or .gpg\n";
print "\n";
print "Options:\n";
print "-h            This help\n";
print "-v            Version number\n";
print "-p password   Set password used to decrypt\n";
print "-d directory  Set directory to recursively decrypt\n";
print "-r            Remove the encrypted files *.pgp or *.gpg\n";
print "-f            Force the use of PGP\n";
print "-g            Force the use of GPG\n";
print "-q            Quietly decrypt\n";
print "\n";
print "Suggested usage:\n\n";
print "> decrypt_data.pl -d directory\n\n";
print "Using GPG 1.2.6\n\n";
print "Enter the password:\n\n";
print "Then, enter the password or decrypting key at the prompt. \n";
print "\n";

}

#
#  'wanted' routine for find
#  Operates on each traversed file
#  If anything ends in .gpg, force usage of gpg
#  Add all .pgp and .gpg files to decrypt list
#
sub visit {
   my $file = $File::Find::name;
   if ( $file =~ /\.gpg$/ ) {
      $gpg = 1;
      push @decrypt_files, $file;
   } elsif ( $file =~ /\.pgp$/ ) {
      push @decrypt_files, $file;
   }
}
