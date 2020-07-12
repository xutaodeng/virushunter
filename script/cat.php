<?php
set_time_limit(0);
function startsWith($haystack, $needle)
{
    return !strncmp($haystack, $needle, strlen($needle));
}

function endsWith($haystack, $needle)
{
    $length = strlen($needle);
    if ($length == 0) {
        return true;
    }

    return (substr($haystack, -$length) === $needle);
}

$stamp = date("Y-m-d-H-i-s");
$q=$_POST["q"];
//echo 'q:'.$q;
$c1 = explode(' ', trim($q));
$path=realpath('..\\tmp');
//realpath('..\\tmp\\export.'.$stamp.'.html');
//echo 'path:'.$path.'<br>';
$exaln = $path.'\\export.'.$stamp.'.html';
$exfa =  $path.'\\export.'.$stamp.'.fa';
$sortaln =  $path.'\\sort.'.$stamp.'.html';
$sortfa =  $path.'\\sort.'.$stamp.'.fa';
$aln = fopen($exaln, "w");
$fa = fopen($exfa, "w");
//echo 'fa path:'.$sortfa.'<br>';

$header = '<html> <head> <link rel="stylesheet" type="text/css" href="../tablestyle.css"></head><script src="../sorttable.js"></script><table class="sortable"><tr><th class="sorttable_nosort">Virus</th><th>Query</th><th>Pair Hit</th><th>Evalue</th><th>LNVNRE</th><th>Identity</th><th class="sorttable_nosort">Alignment</th><th>Query_nt</th></tr>';

fwrite($aln, $header);
fwrite($aln, "\n");
$i=0;
foreach ($c1 as &$value) {
	$i=$i+1;
	if ($i==1){ #evalues
		$e1=$value;
	}
	elseif ($i==2) {
		$e2=$value;
	}
	elseif (trim($value)=='') {
		continue;
	}
	else{
		//echo $value.'<br>';
		$value = str_replace('*', '%2A', $value); 
		#print $value;
		$h1 = fopen('aln\\'.$value.'.html', "r");
		$h2 = fopen('fasta\\'.$value.'.fa', "r");
		while(!feof($h1)){ 
			$buf =fgets($h1);
			if (startsWith($buf,'<tr>') or startsWith($buf,'<td>')){
				fwrite($aln, $buf);
			}
		}
		while(!feof($h2)){ fwrite($fa, fgets($h2));}
		fclose($h1);
		fclose($h2);
	}
}
fwrite($aln, '</table></html>');
fclose($aln);
fclose($fa);
$cmd='E:\\wamp64\\www\\catAlignFA.py '.$e1.' '.$e2.' '.$exaln.' '.$exfa.' '.$sortaln.' '.$sortfa;
// echo $cmd.'<br>';
system($cmd);
echo 'Download Selected <a href="'.'../tmp/sort.'.$stamp.'.html'.'">Hits</a> ';
echo '<a href="'.'../tmp/sort.'.$stamp.'.fa'.'">FASTA</a> ';
?> 
