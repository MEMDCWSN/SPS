#! /usr/bin/perl

$| = 1;

use lib "/home/pbouchard/usr/lib64/perl5/site_perl/5.8.8/x86_64-linux-thread-multi";
use FCGI;
use URI::Escape;

$count = 0;

while (FCGI::accept() >= 0)
{
  $cin = <> . "&" . $ENV{'QUERY_STRING'};

  print("Content-type: text/html\r\n\r\n");

  if (length($cin) > 0)
  {
    @params = split(/&/, $cin);
    foreach $param (@params)
    {
      ($key, $value) = split(/=/, $param);
      $key = uri_unescape($key);
      $value = uri_unescape($value);
      $call .= " " . $key;

      @values = split(/\+/, $value);
      foreach $value (@values)
      {
        $value = quotemeta $value;
        $call .= " " . $value;
      }
    }

    print("<br>\n");
    $cout = `$call 2>&1 1>/dev/null`;
    print("<br>\n");
    print($cout);

    if ($cout ne "")
    {
      sleep 3;
    }

    print("<script language='javascript'>");
    print("  javascript:history.back();");
    print("</script>");
  }
}
