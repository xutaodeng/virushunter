<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01//EN" "http://www.w3.org/TR/html4/strict.dtd">
<html>
<head>
 <meta http-equiv="content-type" content="text/html; charset=ISO-8859-1">
  <title>blast</title>

  <script src="./tree/jquery/jquery.js" type="text/javascript"></script>
  <script src="./tree/jquery/jquery-ui.custom.js" type="text/javascript"></script>
  <script src="./tree/jquery/jquery.cookie.js" type="text/javascript"></script>

  <link href="./tree/src/skin/ui.dynatree.css" rel="stylesheet" type="text/css"> 
  <script src="./tree/src/jquery.dynatree.js" type="text/javascript"></script>
<?php
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
?>
 <script type="text/javascript">
  $(function(){
      $("#tree").dynatree({
      onActivate: function(node) {
	    var db = node.data.title;
	    document.forms["myForm"]["db"].value=db;
        //alert("You activated " + document.forms["myForm"]["db"].value);
      },
      children: [
	{title: "All", isFolder: true, key: "All",
	  children:[
 <?php
	 $handle=opendir(".");
	 $allout=array();
	 $alias=array();
	while ($file = readdir($handle)) 
	{
		if (is_dir($file) && $file!='' && $file!='.' && $file!='..' && $file!='DataTables-1.9.4' && $file!='tree' && $file!='tmp') 
		{
			$children=array();
			$parent="{title:\"NULL\", isFolder: true, children:[";;
			if (is_dir($file.'/blast'))
			{
				if ($handle2=opendir($file.'/blast'))
				{	$c=0;
					while ($file2 = readdir($handle2)) 
					{
						if (endsWith($file2, ".nhr") or endsWith($file2, ".nal")) 
						{
							$fid = substr($file2, 0, strpos($file2, "_blastdb"));
							//{title: "Folder 2", isFolder: true, children: [{title: "Sub-item 2.1"},{title: "Sub-item 2.2"}]},
							if ($fid==$file){
								//echo 'parent? '.$file.'\n';
								$parent="{title:\"".$fid."\", isFolder: true, children:[";
								$rundb = $fid.'\\blast\\'.$fid.'_blastdb';
								array_push($alias,  $rundb);
							}
							else {
								$child='{title: '."\"".$file.'$'.$fid."\"}";
								if (!in_array($child, $children)){
									array_push($children,  $child);
									$c=$c+1;
								}
							}
						}
					}
					if ($c>0){
						$output=$parent.implode(",", $children)."]}\n";
						array_push($allout, $output);
					}
				}
			}
		}
	}
	echo implode(",", $allout);
	#$cmd = "blastdb_aliastool -dblist \"".implode(" ", $alias)."\" -dbtype nucl -out All_blastdb -title \"All_blastdb\""; 
	#echo $cmd;
	#system($cmd, $retval);
	#Undetermined_blastdb Ita-Res-D5-8_blastdb Ita-Res-A6-10_blastdb" -dbtype nucl -out 121116_Phan_Italy_Respiratory_Infection_blastdb -title #"121116_Phan_Italy_Respiratory_Infection_blastdb"
	
 ?>
 	  ]}
      ]
    });
  });
 </script>

<script type="text/javascript">
function clearForms()
{
    document.getElementById('pleaseWait').reset();
}
function validateForm()
{

var x=document.forms["myForm"]["query"].value;
//alert(x);
	var pass=false;
	if (x==null || x=="")
	  {
	  alert("You have empty query sequence!");
	  return false;
	  }
    var radios = document.getElementsByName('blast');
    for (var i = 0; i < radios.length; i++) {
        if (radios[i].checked) {
			pass=true;
			break;
		}
	}
	if (pass==false){
		alert("select your blast program");
		return false;
	}
	// pass=false;
	// var radios2 = document.getElementsByName('db');
    // for (var i = 0; i < radios2.length; i++) {
        // if (radios2[i].checked) {
		// pass=true;
		// break;
		// }
	// }
	// if (pass==false){
	// alert("select your database");
	// return false;
	// }
	document.getElementById('pleaseWait').style.display='block';
	return true;

}
</script>
</head>

  <DIV ID="pleaseWait" STYLE="position:absolute;z-index:5;top:30%;left:35%; display: none">
	<IMG SRC="wait.gif" WIDTH=600 HEIGHT=300>
		  <!-- <IMG SRC="wait.gif" BORDER=1 WIDTH=75 HEIGHT=15><BR><BR>  -->
  </DIV> 
 <body onLoad="clearForms()" onUnload="clearForms()">
  
<form name="myForm" action="blast_run.php" method="post" onsubmit="return validateForm()">

Paste Your Query Sequence: <br>
 <TEXTAREA name="query" rows="10" cols="80"></TEXTAREA><br>
<br>
Select Run to Blast Against: <br>
 <div id="tree"> </div>

<br>
<INPUT type="hidden" name="db" value="All">
<br>
Select blast program: <br>
<input type="radio" name="blast" value="blastn -task megablast">megablast (nucleotide query)<br>
<input type="radio" name="blast" value="blastn -task blastn">blastn (nucleotide query)<br>
<input type="radio" name="blast" value="tblastx">tblastx (nucleotide query)<br>
<input type="radio" name="blast" value="tblastn">tblastn (protein query)<br>
Evalue: <input type="text" name="evalue" value="10"><br>

<input type="submit" name = "submit" value="blast" style="background-color:lightgreen">

</form> 
</body>
</html>