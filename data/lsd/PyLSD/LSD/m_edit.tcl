#    file: m_edit.tcl
#    Copyright (C) 2000 CNRS, UMR 7312, Jean-Marc Nuzillard.
#
#    This file is part of LSD.
#
#    LSD is free software; you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation; either version 2 of the License, or
#    (at your option) any later version.
#
#    LSD is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with LSD; if not, write to the Free Software
#    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

if {$argc == 0} {
	set fileName Untitled
} elseif {$argc == 1} {
	set fileName [lindex $argv 0]
} else {
	puts {Usage: wish m_edit [filename]}
	exit 1
}

set fileType {
	{"LSD Coordinates"	{.coo}	}
	{"SDF files"		{.sdf}	}
}

set Radius 12
set Shift 3
set curElem C
set curCharge 0

set cWidth 700
set cHeight 500
set lineWidth 2
set fontSize 12
set numFont [ font create -size $fontSize -weight bold -family Helvetica ]
set molTitleFont [ font create -size [ expr $fontSize * 3 ] -weight bold -family Helvetica ]
set xTitle [ expr $cWidth * 0.9 ]
set yTitle [ expr $cHeight * 0.1 ]

frame .mbar -relief raised -bd 2
canvas .c -width $cWidth -height $cHeight -bg white
pack .mbar .c -side top -fill x
menubutton .mbar.file -text File -menu .mbar.file.menu
menubutton .mbar.view -text View -menu .mbar.view.menu
menubutton .mbar.select -text Select -menu .mbar.select.menu
menubutton .mbar.buffer -text Buffer -menu .mbar.buffer.menu
menubutton .mbar.arrange -text Arrange -menu .mbar.arrange.menu
pack .mbar.file .mbar.view .mbar.select .mbar.buffer .mbar.arrange -side left
menu .mbar.file.menu
# .mbar.file.menu add command -label "New" -command "newFile"
.mbar.file.menu add command -label "Open" -command "openFile"
.mbar.file.menu add command -label "Save" -command "saveFile"
.mbar.file.menu add command -label "Save As" -command "saveAsFile"
.mbar.file.menu add command -label "Exit" -command "exitFile"
menu .mbar.view.menu
.mbar.view.menu add command -label "Next" -command "nextMol"
.mbar.view.menu add command -label "Previous" -command "prevMol"
# .mbar.view.menu add command -label "New" -command "newMol"
menu .mbar.select.menu
.mbar.select.menu add command -label "None" -command "unSelectAll"
.mbar.select.menu add command -label "All" -command "selectAll"
.mbar.select.menu add command -label "Invert Current" -command "invCurSel"
.mbar.select.menu add command -label "Invert All" -command "invAllSel"
.mbar.select.menu add command -label "Keep" -command "keepSelectedOnly"
menu .mbar.buffer.menu
.mbar.buffer.menu add command -label "Save to" -command "saveToBuffer"
.mbar.buffer.menu add command -label "Load from" -command "loadFromBuffer"
menu .mbar.arrange.menu
.mbar.arrange.menu add command -label "Horiz. Flip" -command "hFlip"
.mbar.arrange.menu add command -label "Vert. Flip" -command "vFlip"

proc mkAtom {i x y elt label charge} {
	global atX atY atElt atLabel atCharge atId atBond idAt numFont

	set atX($i) $x
	set atY($i) $y
	set atElt($i) $elt
	set atLabel($i) $label
	set atCharge($i) $charge
	if {$charge > 0} {
		if {$charge == 1} {
			set txtCharge " +"
		} else {
			set txtCharge " $charge+"
		}
	} elseif {$charge < 0} {
		if {$charge == -1} {
			set txtCharge " -"
		} else {
			set txtCharge " [expr {abs($charge)}]-"
		}
	} else {
		set txtCharge ""
	}
	set txt $elt$label$txtCharge
	set atBond($i) {}
	set id [.c create text $x $y -text $txt -font $numFont -tags node]
	set atId($i) $id
	set idAt($id) $i
}

proc hFlip {} {
	global atX atElt cWidth modified

	foreach i [array names atElt] {
		set atX($i) [ expr $cWidth - $atX($i) ]
	}
	set modified 1
	saveMol
	drawMol
}

proc vFlip {} {
	global atY atElt cHeight modified

	foreach i [array names atElt] {
		set atY($i) [ expr $cHeight - $atY($i) ]
	}
	set modified 1
	saveMol
	drawMol
}

proc atToBond {a1 a2} {
	expr {($a1 < $a2) ? "$a1,$a2" : "$a2,$a1"}
}

proc defBond order {
	global firstAtom idAt modified

	set curAtom [.c find withtag current]
	if {($firstAtom != "") && ($curAtom != "")} {
		set a1 $idAt($firstAtom)
		set a2 $idAt($curAtom)
		mkBond $a1 $a2 $order
	}
	set modified 1
}

proc mkBond {a1 a2 order} {
	global atX atY atBond bondOrder bondId idBond lineWidth

	set bond [atToBond $a1 $a2]
	rmAllBonds $bond
	lappend atBond($a1) $a2
	lappend atBond($a2) $a1
	set x1 $atX($a1)
	set y1 $atY($a1)
	set x2 $atX($a2)
	set y2 $atY($a2)
	set cooList [mkCooList $x1 $y1 $x2 $y2 $order]
	set newList {}
	for {set i 0} {$i < $order} {incr i} {
		set coo [lindex $cooList $i]
		set new [eval .c create line $coo -width $lineWidth -tags bond]
		set idBond($new) $bond
		.c lower $new
		lappend newList $new
	}
	set bondOrder($bond) $order
	set bondId($bond) $newList
}

proc mkCooList {x1 y1 x2 y2 order} {
	switch $order {
		1 	{
			set c1 [mkCoo $x1 $y1 $x2 $y2 0]
			return [list "$c1"]
		}
		2	{
			set c1 [mkCoo $x1 $y1 $x2 $y2 0.5]
			set c2 [mkCoo $x1 $y1 $x2 $y2 -0.5]
			return [list "$c1" "$c2"]
		}
		3	{
			set c1 [mkCoo $x1 $y1 $x2 $y2 1]
			set c2 [mkCoo $x1 $y1 $x2 $y2 0]
			set c3 [mkCoo $x1 $y1 $x2 $y2 -1]
			return [list "$c1" "$c2" "$c3"]
		}
	}
}

proc mkCoo {x1 y1 x2 y2 sp} {
	global Radius Shift

	set split [expr $Shift * $sp]
	set dx [expr $x2 - $x1]
	set dy [expr $y2 - $y1]
	set d [expr sqrt($dx * $dx + $dy * $dy)]
	set ux [expr $dx / $d]
	set uy [expr $dy / $d]
	set x3 [expr $x1 + $Radius * $ux - $split * $uy]
	set y3 [expr $y1 + $Radius * $uy + $split * $ux]
	set x4 [expr $x2 - $Radius * $ux - $split * $uy]
	set y4 [expr $y2 - $Radius * $uy + $split * $ux]
	return "$x3 $y3 $x4 $y4"
}

proc mkTitle {} {
	global xTitle yTitle molTitle molTitleFont molSelected

	.c delete title
	set color [ expr { $molSelected ? "black" : "red" } ]
	.c create text $xTitle $yTitle -text $molTitle -font $molTitleFont -tags title -fill $color
}

proc unSelectAll {} {
	setSelectionToAll 0
}

proc selectAll {} {
	setSelectionToAll 1
}

proc setSelectionToAll v {
	global mol molSelected numMol

	for {set i 1} {$i <= $numMol} {incr i} {
		set mol($i) [ lreplace $mol($i) 5 5 $v ]
	}
	set molSelected $v
	mkTitle
}

proc invCurSel {} {
	global molSelected

	set molSelected [ expr 1 - $molSelected ]
	mkTitle
}

proc invAllSel {} {
	global mol molSelected numMol

	for {set i 1} {$i <= $numMol} {incr i} {
		set mol($i) [ lreplace $mol($i) 5 5 [ expr 1 - [ lindex $mol($i) 5 ] ] ]
	}
	invCurSel
}

proc keepSelectedOnly {} {
	global mol numMol curMol modified

	saveMol

	set numKept 0
	set numRemoved 0
	set curMol2 $curMol

	for {set i 1} {$i <= $numMol} {incr i} {
		set selected [ lindex $mol($i) 5 ]
		if { $selected } {
			incr numKept
			set mol2($numKept) $mol($i)
		} {
			incr numRemoved
			if { $i <= $curMol } {
				incr curMol2 -1
			}
		}
	}
	if { [ expr $numRemoved && [ okToRemove $numRemoved $numKept ] ] } {
		clearMol
		set modified 1
		set numMol $numKept
		set curMol [ expr { $curMol2 ? $curMol2 : 1 } ]
		if { $numMol } {
			unset mol
			for {set i 1} {$i <= $numMol} {incr i} {
				set mol($i) $mol2($i)
			}
			drawMol
		} {
			set curMol 1
			set numMol 1
		}
		showCurMol
	}
}

proc okToRemove { nrem nkept } {
	return [tk_dialog .okrem {Keep selected structures only} \
		"Keeping $nkept structure(s)\nand removing $nrem structure(s). OK ?" \
		warning 0 {No} {Yes} ]
}

proc saveToBuffer {} {
	global buffer mol

	saveMol
	catch {unset buffer}
	foreach i [array names mol] {
		set buffer($i) $mol($i)
	}
}

proc loadFromBuffer {} {
	global buffer mol modified curMol numMol

	if { [info exists buffer] } {
		set ok 1
		if { [ expr $modified && [tk_dialog .d {Proceed} "Proceed ?" \
		warning 1 {Yes} {No} ] ] } {
			set ok 0
		}
		if { $ok } {
			catch {unset mol}
			foreach i [array names buffer] {
				set mol($i) $buffer($i)
			}
			set numMol [array size mol]
			set curMol 1
			showCurMol
			drawMol
			set modified 0
		}
	}
}

proc moveAtom {id xDist yDist} {
	global idAt atX atY atBond bondOrder bondId

	.c move $id $xDist $yDist
	set a1 $idAt($id)
	incr atX($a1) $xDist
	incr atY($a1) $yDist
	set x1 $atX($a1)
	set y1 $atY($a1)
	foreach a2 $atBond($a1) {
		set bond [atToBond $a1 $a2]
		set order $bondOrder($bond)
		set lignes $bondId($bond)
		set x2 $atX($a2)
		set y2 $atY($a2)
		set coos [mkCooList $x1 $y1 $x2 $y2 $order]
		for {set i 0} {$i < $order} {incr i} {
			set coo [lindex $coos $i]
			set ligne [lindex $lignes $i]
			eval .c coords $ligne $coo
		}
	}
}

proc saveMol {} {
	global atX atY atElt atLabel atCharge atId atBond idAt
	global bondOrder bondId curMol molTitle molSelected mol

	set numAtom [array size atElt]
	set numBond [array size bondOrder]
	set listAtom ""
	foreach i [array names atElt] {
		lappend listAtom [list $i $atX($i) $atY($i) $atElt($i) $atLabel($i) $atCharge($i)]
	}
	set listBond ""
	foreach i [array names bondOrder] {
		lappend listBond [list $i $bondOrder($i)]
	}
	set mol($curMol) [list $numAtom $numBond $listAtom $listBond $molTitle $molSelected]
}

proc clearMol {} {
	global atX atY atElt atLabel atCharge atId atBond idAt
	global bondOrder bondId idBond molTitle molSelected nextAtom

	.c delete node
	.c delete bond
	.c delete title
	catch {unset atX atY atElt atLabel atCharge atId atBond idAt bondOrder bondId idBond molTitle molSelected}
	set nextAtom 0
}

proc drawMol {} {
	global mol nextAtom
	global curMol molTitle molSelected

	clearMol
	set molInfo $mol($curMol)
	set numAtom [lindex $molInfo 0]
	set numBond [lindex $molInfo 1]
	set atomList [lindex $molInfo 2]
	set bondList [lindex $molInfo 3]
	set molTitle [lindex $molInfo 4]
	set molSelected [lindex $molInfo 5]
	set nextAtom $numAtom
	foreach atom $atomList {
		eval mkAtom $atom
	}
	foreach bond $bondList {
		set str [lindex $bond 0]
		set order [lindex $bond 1]
		set lis [split $str , ]
		set a1 [lindex $lis 0]
		set a2 [lindex $lis 1]
		mkBond $a1 $a2 $order
	}
	mkTitle
}

proc newMol {} {
	global curMol numMol modified

	saveMol
	clearMol
	incr numMol
	set curMol $numMol
	showCurMol
	set modified 1
}

proc showCurMol {} {
	global fileName curMol numMol

	set text "m_edit - [file tail $fileName] $curMol / $numMol"
	wm title . $text
}

proc nextMol {} {
	global numMol curMol

	if {$numMol == $curMol} return
	saveMol
	clearMol
	incr curMol
	showCurMol
	drawMol
}

proc prevMol {} {
	global numMol curMol

	if {$curMol == 1} return
	saveMol
	clearMol
	incr curMol -1
	showCurMol
	drawMol
}

proc rmAtom id {
	global idAt atX atY atElt atLabel atCharge atId idAt atBond

	.c delete $id
	set a1 $idAt($id)
	unset atX($a1) atY($a1) atElt($a1) atLabel($a1) atCharge($a1)
	unset atId($a1) idAt($id)
	if {[info exists atBond($a1)]} {
		foreach a2 $atBond($a1) {rmAllBonds [atToBond $a1 $a2]}
		unset atBond($a1)
	}
}

proc rmOneBond bond {
	global bondOrder

	set order [expr $bondOrder($bond) - 1]
	rmAllBonds $bond
	if {$order != 0} {
		set atoms [split $bond ,]
		set a1 [lindex $atoms 0]
		set a2 [lindex $atoms 1]
		mkBond $a1 $a2 $order
	}
}

proc rmAllBonds bond {
	global bondId bondOrder atBond

	if {[info exists bondId($bond)]} {
		set atoms [split $bond ,]
		set a1 [lindex $atoms 0]
		set a2 [lindex $atoms 1]
		foreach id $bondId($bond) {.c delete $id}
		unset bondId($bond) bondOrder($bond)
		listRmOne $atBond($a1) $a2
		listRmOne $atBond($a2) $a1
	}
}

proc listRmOne {list el} {
	set i [lsearch -exact $list $el]
	if {$i >= 0} {
		lreplace $list $i $i
	} else {
		set $list
	}
}

proc centerCoo {xtab ytab n} {
	upvar $xtab atx
	upvar $ytab aty

	set xmin $atx(0)
	set ymin $aty(0)
	set xmax $xmin
	set ymax $ymin
	for {set a 1} {$a < $n} {incr a} {
		set x $atx($a)
		set y $aty($a)
		if {$x > $xmax} {set xmax $x}
		if {$y > $ymax} {set ymax $y}
		if {$x < $xmin} {set xmin $x}
		if {$y < $ymin} {set ymin $y}
	}
	set xmean [expr ($xmin + $xmax)/2]
	set ymean [expr ($ymin + $ymax)/2]
	for {set a 0} {$a < $n} {incr a} {
		set atx($a) [expr $atx($a) - $xmean]
		set aty($a) [expr $aty($a) - $ymean]
	}
	set dx [expr $xmax - $xmin]
	set dy [expr $ymax - $ymin]
	list $dx $dy
}

proc scaleCenterDisplay numAtom {
	global atX atY cWidth cHeight

	set delta [centerCoo atX atY $numAtom]
	set dx [lindex $delta 0]
	set dy [lindex $delta 1]
	set cntrX [expr $cWidth/2]
	set cntrY [expr $cHeight/2]
	set echx [expr 0.8 * $cWidth / $dx]
	set echy [expr 0.8 * $cHeight / $dy]
	set ech [expr {($echx < $echy) ? $echx : $echy}]
	for {set i 0} {$i < $numAtom} {incr i} {
		set x [expr $cntrX+$atX($i)*$ech]
		set y [expr $cntrY-$atY($i)*$ech]
		set atX($i) [expr round($x)]
		set atY($i) [expr round($y)]
	}
}

proc openFile {} {
	global fileName fileType

	set fn [tk_getOpenFile -filetypes $fileType -parent .c]
	if [string compare $fn ""] {
		closeFile
		set fileName $fn
		readFile
	}
}

proc readFile {} {
	global fileName modified curMol

	set ext [file extension $fileName]
	switch $ext {
		.coo {readFileCoo}
		.sdf {readFileMol}
	}
	set curMol 1
	showCurMol
	drawMol
	set modified 0
}

proc readFileCoo {} {
	global fileName
	global atX atY atElt atLabel atCharge bondOrder mol
	global numMol curMol molTitle molSelected

	set curMol 0
	set f [open $fileName r]
	gets $f line
	while {[gets $f line] >= 0} {
		set numAtom [lindex $line 0]
		set molTitle [lindex $line 1]
		incr curMol
		for {set a1 0} {$a1 < $numAtom} {incr a1} {
			gets $f line
			set atX($a1) [lindex $line 17]
			set atY($a1) [lindex $line 18]
			set atElt($a1) [lindex $line 2]
			set atCharge($a1) [lindex $line 3]
			set atLabel($a1)  [lindex $line 1]
			set nbond  [lindex $line 4]
			for {set i 0} {$i < $nbond} {incr i} {
				set j [expr 5+2*$i]
				set k [expr $j + 1]
				set a2 [lindex $line $j]
				if {$a1 < $a2} {
					set bond "$a1,$a2"
					set order [lindex $line $k]
					set bondOrder($bond) $order
				}
			}
		}
		scaleCenterDisplay $numAtom
		set molSelected 1
		saveMol
		unset atX atY atElt atLabel atCharge bondOrder
	}
	set numMol $curMol
}

proc readFileMol {} {
	global fileName
	global atX atY atElt atLabel atCharge bondOrder mol cHeight cWidth
	global numMol curMol molTitle molSelected

	set curMol 0
	set f [open $fileName r]
	set status 1
	while {$status >= 0} {
		incr curMol
		set molTitle $curMol
		gets $f line
		gets $f line
		gets $f line
		if {$curMol == 1} {gets $f line}
		set numAtom [lindex $line 0]
		set numBond [lindex $line 1]
		set numEndline [lindex $line 10]
		for {set a 0} {$a < $numAtom} {incr a} {
			gets $f line
			set atX($a) [lindex $line 0]
			set atY($a) [lindex $line 1]
			set atElt($a) [lindex $line 3]
			set c [lindex $line 5]
			if {$c == 0} {
				set atCharge($a) 0
			} else {
				set atCharge($a) [expr 4 - $c]
			}
			set atLabel($a) [expr $a + 1]
		}
		for {set i 0} {$i < $numBond} {incr i} {
			gets $f line
			set a1 [expr [lindex $line 0] - 1]
			set a2 [expr [lindex $line 1] - 1]
			set bondOrder([atToBond $a1 $a2]) [lindex $line 2]
		}
		for {set i 0} {$i < $numEndline} {incr i} {
			gets $f line
		}
		gets $f line
		scaleCenterDisplay $numAtom
		set molSelected 1
		saveMol
		unset atX atY atElt atLabel atCharge bondOrder
		set status [gets $f line]
	}
	set numMol $curMol
}

proc saveFile {} {
	global fileName

	if [string compare "$fileName" "Untitled"] {
		writeFile
	} else {
		set choice [tk_dialog .dd {Improper filename} \
{Cannot save "Untitled" file} warning 1 {Save As} {Cancel}]
		if {!$choice} {
			saveAsFile
		}
	}
}

proc saveAsFile {} {
	global fileName fileType

	set fn [tk_getSaveFile -filetypes $fileType -parent .c \
-initialfile Untitled -defaultextension .coo]
	if [string compare $fn ""] {
		set fileName $fn
		saveFile
		showCurMol
	}
}

proc writeFile {} {
	global fileName modified

	set ext [file extension $fileName]
	switch $ext {
		.coo {writeFileCoo}
		.sdf {writeFileMol}
	}
	set modified 0
}

proc writeFileCoo {} {
	global fileName
	global mol curMol numMol cHeight
	
	saveMol
	set f [open $fileName w]
	puts $f DRAW
	for {set i 1} {$i <= $numMol} {incr i} {
		set molInfo $mol($i)
		set numAtom [lindex $molInfo 0]
		set atomList [lindex $molInfo 2]
		set bondList [lindex $molInfo 3]
		set molTitle [lindex $molInfo 4]
		puts $f [format "%3d%5d" $numAtom $molTitle]
		foreach atomInfo $atomList {
			set index [lindex $atomInfo 0]
			set atx($index) [lindex $atomInfo 1]
			set aty($index) [expr $cHeight - [lindex $atomInfo 2]]
			set atelt($index) [lindex $atomInfo 3]
			set atlabel($index) [lindex $atomInfo 4]
			set atcharge($index) [lindex $atomInfo 5]
			set numnei($index) 0
			for {set j 1} {$j <= 6} {incr j} {
				set nei($index,$j) 0
				set ord($index,$j) 0
			}
		}
		foreach bondInfo $bondList {
			set bond [lindex $bondInfo 0]
			set order [lindex $bondInfo 1]
			set atoms [split $bond ,]
			set a1 [lindex $atoms 0]
			set a2 [lindex $atoms 1]
			set ia1 [incr numnei($a1)]
			set nei($a1,$ia1) $a2
			set ord($a1,$ia1) $order
			set ia2 [incr numnei($a2)]
			set nei($a2,$ia2) $a1
			set ord($a2,$ia2) $order
		}
		centerCoo atx aty $numAtom
		for {set a 0} {$a < $numAtom} {incr a} {
			set line1 [format "%3d%3d%3s%3d%2d" \
$a $atlabel($a) $atelt($a) $atcharge($a) $numnei($a)]
			set line2 "  "
			for {set j 1} {$j <= 6} {incr j} {
				set line2 $line2[format "%3d%2d" \
$nei($a,$j) $ord($a,$j)]
			}
			set line3 [format "%10d %10d" \
$atx($a) $aty($a)]
			puts $f "$line1$line2$line3"
		}
	}
	flush $f
}

proc writeFileMol {} {
	global fileName
	global mol curMol numMol cHeight
	
	saveMol
	set f [open $fileName w]
	for {set i 1} {$i <= $numMol} {incr i} {
		set numCharge 0
		puts $f [format "\n\n"]
		set molInfo $mol($i)
		set numAtom [lindex $molInfo 0]
		set numBond [lindex $molInfo 1]
		set atomList [lindex $molInfo 2]
		set bondList [lindex $molInfo 3]
		foreach atomInfo $atomList {
			set index [lindex $atomInfo 0]
			set atx($index) [lindex $atomInfo 1]
			set aty($index) [expr $cHeight - [lindex $atomInfo 2]]
			set atelt($index) [lindex $atomInfo 3]
			set atcharge($index) [lindex $atomInfo 5]
		}
		for {set a 0} {$a < $numAtom} {incr a} {
			if {$atcharge($a) != 0} {
				incr numCharge
			}
		}
#		puts $f [format "%3d%3d  0  0  0  0  0  0  0  0%3d V2000" \
#$numAtom $numBond [expr {$numCharge ? 2 : 1}]]
		puts $f [format "%3d%3d  0  0  0  0  0  0  0  0999 V2000" $numAtom $numBond]
		centerCoo atx aty $numAtom
		set dsum 0
		foreach bondInfo $bondList {
			set bond [lindex $bondInfo 0]
			set order [lindex $bondInfo 1]
			set atoms [split $bond ,]
			set a1 [lindex $atoms 0]
			set a2 [lindex $atoms 1]
			set x [expr $atx($a2) - $atx($a1)]
			set y [expr $aty($a2) - $aty($a1)]
			set dsum [expr hypot($x, $y) + $dsum]
		}
		set ech [expr 1.0 / ($dsum/$numBond)]
		for {set a 0} {$a < $numAtom} {incr a} {
			set x [expr $atx($a) * $ech]
			set y [expr $aty($a) * $ech]
			puts $f [format "%10.4f%10.4f%10.4f %-3s\
  0%3d  0  0  0  0  0  0  0  0  0  0" $x $y 0.0 $atelt($a) [expr {$atcharge($a) ? [expr {4 - $atcharge($a)}] : 0}]]
		}
		foreach bondInfo $bondList {
			set bond [lindex $bondInfo 0]
			set order [lindex $bondInfo 1]
			set atoms [split $bond ,]
			set a1 [expr [lindex $atoms 0] + 1]
			set a2 [expr [lindex $atoms 1] + 1]
			puts $f [format "%3d%3d%3d  0  0  0  0" $a1 $a2 $order]
		}
		if {$numCharge} {
			set part1 [format "M  CHG%3d" $numCharge]
			set part2 ""
			for {set a 0} {$a < $numAtom} {incr a} {
				if {$atcharge($a) != 0} {
					set part2 $part2[format "%4d%4d" [expr $a + 1] $atcharge($a)]
				}
			}
			puts $f "$part1$part2"
		}
		puts $f [format "M  END"]
		puts $f \$\$\$\$
	}
	flush $f
}

proc closeFile {} {
	global modified mol fileName numMol curMol nextAtom

	set cancelled 0
	if {$modified} {
		set choice [tk_dialog .d {File modified} \
"File \"$fileName\" was modified. Do you want to save it?" \
warning 3 {Save} {Save As} {Discard} {Cancel}]
		switch $choice {
		0	{saveFile}
		1	{saveAsFile}
		2	{}
		3	{set cancelled 1}
		}
	}
	if {!$cancelled} {
		catch {unset mol}
		set fileName Untitled
		set numMol 1
		set curMol 1
		set nextAtom 0
		set modified 0
		showCurMol
		clearMol
	}
	set cancelled
}

proc newFile {} {
	closeFile
}

proc exitFile {} {
	if {![closeFile]} {exit 0}
}

bind .c <Button-3> {
	global nextAtom curElem curCharge

	mkAtom $nextAtom %x %y $curElem [incr nextAtom] $curCharge
	set modified 1
}
bind .c 0 {
	global firsAtom

	set firstAtom [.c find withtag current]
}
bind .c 1 {defBond 1}
bind .c 2 {defBond 2}
bind .c 3 {defBond 3}

bind .c <BackSpace> {
	global idBond idAt

	set curObj [.c find withtag current]
	if {[info exists idAt($curObj)]} {rmAtom $curObj}
	if {[info exists idBond($curObj)]} {
		eval rmOneBond $idBond($curObj)
	}
	set modified 1
}

.c bind node <Any-Enter> {
	.c itemconfigure current -fill red
}
.c bind node <Any-Leave> {
	.c itemconfigure current -fill black
}
.c bind bond <Any-Enter> {
	.c itemconfigure current -fill red
}
.c bind bond <Any-Leave> {
	.c itemconfigure current -fill black
}
.c bind node <Button-1> {
	set curX %x
	set curY %y
}
.c bind node <B1-Motion> {
	set dx [expr %x-$curX]
	set dy [expr %y-$curY]
	moveAtom [.c find withtag current] $dx $dy
	set curX %x
	set curY %y
	set modified 1
}
.c bind title <Button-1> {
	invCurSel
}

set numMol 1
set curMol 1
set nextAtom 0
set modified 0

if {$argc == 1} {
	readFile
}
showCurMol

focus .c
