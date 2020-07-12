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
$c1 = explode(' ', trim($q));
#echo $q;
$i=0;
$rval='';
foreach ($c1 as &$value) {
	$i=$i+1;
	if ($i==1){ #evalues
		$e1=floatval($value);
	}
	elseif ($i==2) {
		$e2=floatval($value);
	}
	elseif (trim($value)=='') {
		continue;
	}
	else{
		$value = str_replace('*', '%2A', $value); 
		$h1 = fopen('aln\\'.$value.'.html', "r");
		$h2 = fopen('fasta\\'.$value.'.fa', "r");
		$w1 = fopen('aln\\'.$value.'.tmp.html', "w");
		$w2 = fopen('fasta\\'.$value.'.tmp.fa', "w");
		$count=0;
		while(!feof($h1)){ 
			$buf =fgets($h1);
			if (startsWith($buf,'<td>') or startsWith($buf,'<tr>')){
				$arr=explode(' ', str_replace('<td>', ' ', $buf));
				$e = floatval($arr[4]);
				if ($e1<=$e and $e<=$e2)
				{
					$fa_header =fgets($h2);
					$fa_seq =fgets($h2);
					fwrite($w1, $buf);
					fwrite($w2, $fa_header);
					fwrite($w2, $fa_seq);
					$count=$count+1;
				}
			}
			else{
				fwrite($w1, $buf);
			}
		}
		fclose($h1);
		fclose($h2);
		fclose($w1);
		fclose($w2);
		$rval=$rval.(string) $count.' ';
	}
}

echo $rval;
?> 
