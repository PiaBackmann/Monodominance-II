extensions [table]

globals [
  worldarea
  maxdispersal
  seedpool
  specieslist
  maxpft

  no_died
  tab-dispersalrange
  tab-seednumber-perpatch
  tab-speciescolor
  tab-startseeds
  tab-shadetolerance
  tab-similaradjustment
  tab-sameadjustment

 ; disturbance-range
  disturbance-severity
  mortality-rate

  clusteroutputlist
]




;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Setup Procedures ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

breed [trees tree]

trees-own[
  species
  speciesnumber
  status
  dispersalrange
  dispersalarea
  seednumberperpatch
  shadetolerance
  similaradjustment
  sameadjustment
]

patches-own[
  former_species
  former_shadetolerance
]


to setup
  __clear-all-and-reset-ticks
  resize-world 0 (worldsize - 1) 0 (worldsize - 1)
  set-patch-size round(450 / worldsize)
  set worldarea count patches
  set-default-shape trees "square"
  ask patches [ set pcolor brown - 2 ]
  ask patches [ set former_shadetolerance 50]
  set-parameter
  calculate-maxdispersal
  seed-randomly
  create-plot
end



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Parameter Setting ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to set-parameter

  set worldsize WorldSize
  set disturbance-severity 30.000000
  set mortality-rate 0.015000

 ; set clusteroutputlist [1 250 500 750 1000 1500 2000 2500 3000 3500 4000 4500 5000 5500 6000 6500 7000 7500 8000 8500 9000 9500 10001]

  ;;worldsize can be found below in the code for the bar. Should be changed.

  set tab-dispersalrange table:make
  table:put tab-dispersalrange "pft1" 1
  table:put tab-dispersalrange "pft2" 20
  table:put tab-dispersalrange "pft3" 20
  table:put tab-dispersalrange "pft4" 20
  table:put tab-dispersalrange "pft5" 20
  table:put tab-dispersalrange "pft6" 20
  table:put tab-dispersalrange "pft7" 20
  table:put tab-dispersalrange "pft8" 20

  set tab-seednumber-perpatch table:make
  table:put tab-seednumber-perpatch "pft1"  seedmass
  table:put tab-seednumber-perpatch "pft2"  20
  table:put tab-seednumber-perpatch "pft3"  20
  table:put tab-seednumber-perpatch "pft4"  20
  table:put tab-seednumber-perpatch "pft5"  20
  table:put tab-seednumber-perpatch "pft6"  20
  table:put tab-seednumber-perpatch "pft7"  20
  table:put tab-seednumber-perpatch "pft8"  20

  set tab-startseeds table:make          ;;Summe sollte 100 sein
  table:put tab-startseeds "pft1" 12.5
  table:put tab-startseeds "pft2" 12.5
  table:put tab-startseeds "pft3" 12.5
  table:put tab-startseeds "pft4" 12.5
  table:put tab-startseeds "pft5" 12.5
  table:put tab-startseeds "pft6" 12.5
  table:put tab-startseeds "pft7" 12.5
  table:put tab-startseeds "pft8" 12.5

  set tab-speciescolor table:make
  table:put tab-speciescolor "pft1" 0    ;; 0   = "black"
  table:put tab-speciescolor "pft2" 102    ;; 63  = "green"
  table:put tab-speciescolor "pft3" 103    ;; 105 = "blue"
  table:put tab-speciescolor "pft4" 104    ;; 16  = "red"
  table:put tab-speciescolor "pft5" 105    ;; 35   = "brown"
  table:put tab-speciescolor "pft6" 106   ;; 63  = "green"
  table:put tab-speciescolor "pft7" 107    ;; 105 = "blue"
  table:put tab-speciescolor "pft8" 109   ;; 16  = "red"

  ;; Describes how shade-tolerant different species are supposed to be
  ;; Influences the amount of adjustment through "similaradjustment"
  ;; Between 1 and 100
  set tab-shadetolerance table:make
  table:put tab-shadetolerance "pft1" 100
  table:put tab-shadetolerance "pft2"  50
  table:put tab-shadetolerance "pft3"  50
  table:put tab-shadetolerance "pft4"  50
  table:put tab-shadetolerance "pft5" 50
  table:put tab-shadetolerance "pft6"  50
  table:put tab-shadetolerance "pft7"  50
  table:put tab-shadetolerance "pft8"  50

  ;; How strongly does the species formerly inhabiting the plot
  ;; influence the chances of a similar species to inhabit that plot?
  ;; Can be positive or negative. default is 0.
  set tab-similaradjustment table:make
  table:put tab-similaradjustment "pft1" 0.000000
  table:put tab-similaradjustment "pft2" 0.000000
  table:put tab-similaradjustment "pft3" 0.000000
  table:put tab-similaradjustment "pft4" 0.000000
  table:put tab-similaradjustment "pft5" 0.000000
  table:put tab-similaradjustment "pft6" 0.000000
  table:put tab-similaradjustment "pft7" 0.000000
  table:put tab-similaradjustment "pft8" 0.000000

  ;; How strongly does the species formerly inhabiting the plot
  ;; influence the chances of the same species to inhabit that plot?
  ;; Can be positive or negative. default is 0.
  set tab-sameadjustment table:make
  table:put tab-sameadjustment "pft1" 0.000000
  table:put tab-sameadjustment "pft2" 0.000000
  table:put tab-sameadjustment "pft3" 0.000000
  table:put tab-sameadjustment "pft4" 0.000000
  table:put tab-sameadjustment "pft5" 0.000000
  table:put tab-sameadjustment "pft6" 0.000000
  table:put tab-sameadjustment "pft7" 0.000000
  table:put tab-sameadjustment "pft8" 0.000000

  set specieslist (table:keys tab-seednumber-perpatch)
  set maxpft (length specieslist)

end

to calculate-maxdispersal
  ;;maxima aller dispersal-ranges bestimmen
  set maxdispersal 0
  foreach ( table:keys tab-dispersalrange )
  [
    [x] -> let drange table:get tab-dispersalrange x
    if (drange > maxdispersal) [set maxdispersal drange]
  ]
end

to seed-randomly
   foreach (table:keys tab-startseeds)
   [
    [x] -> ask n-of floor ((table:get tab-startseeds x * worldarea) / 100) patches with [not any? trees-here]
       [ sprout-trees 1 [set-species x] ]
   ]

end

to set-species [species-name]
  set species species-name
  set speciesnumber (position species-name specieslist) + 1
  set color table:get tab-speciescolor species-name
  set dispersalrange table:get tab-dispersalrange species-name
  set seednumberperpatch table:get tab-seednumber-perpatch species-name
  set dispersalarea (count patches in-radius dispersalrange)
  set shadetolerance table:get tab-shadetolerance species-name
  set similaradjustment table:get tab-similaradjustment species-name
  set sameadjustment table:get tab-sameadjustment species-name
  set status "juvenile"
end


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Runtime Procedures ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


to go

   if (ticks = 10000) [export-world ("arena.out") stop]

   ask patches [ set pcolor brown - 2 ]
   set no_died 0
   mortality
   disturbance
   recruitment
   update-status

;   foreach clusteroutputlist
;  [
;    [x] -> if ticks = x
;           [
;                export-view (word "view_" behaviorspace-run-number "_" x"_.png")
;                export-world (word "world_" behaviorspace-run-number "_" x"_.out")
;           ]
;   ]
   tick
end


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; mortality ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
to mortality
  ask trees [

    if (random 100 < (mortality-rate * 100) )
    [ ;;Baum stirbt
      ask patch-here [ set pcolor 23 ]
      ;ask patch-here [ set pcolor brown - 2 ]
      set status "dead"
    ]
  ]
end


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; recruitment ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to recruitment
  ask patches with [not any? trees-here] ;; Schleife über alle leeren Patches
 ; ask patches with [pcolor = 23 ]
  [
    set seedpool []
    foreach (table:keys tab-seednumber-perpatch) [ set seedpool lput 0 seedpool ]  ;;seedpool mit 0 initialisieren

    ask trees in-radius maxdispersal
    [
      if (distance myself <= dispersalrange and status = "alive")
      [
        let localseeds seednumberperpatch
        ;;Prozess 1: ähnliche Schattentoleranz bevorzugen
        ;;  set localseeds localseeds * ( 1 + similaradjustment * ( 1 - abs(shadetolerance - former_shadetolerance) / 100 ) )
        ;;Prozess 2: gleiche Art bevorzugen
       ;;   if ( species = former_species) [ set localseeds localseeds * (1 + sameadjustment) ]

        ;; Take shade tolerance and shading into account

        ;; Take Janzen Connell into account

        ifelse (localseeds >= 0)
        [
          let index (speciesnumber - 1)
          let new-value ( (item index  seedpool) + localseeds )
          set seedpool replace-item index seedpool new-value
        ]
        [
          show "Fehler!!!"
        ]

      ]
    ] ;;Ende Baum-Schleife

    ;;neue Art auswürfeln
    if (sum seedpool > 0)
    [
      let rand (random sum seedpool)
      let id 0
      let rest-seedpool 0
      let new-species-no 0
      let runloop True
      while [runloop]
      [
        if (rand < ( (item id seedpool) + rest-seedpool ) )
        [
          set new-species-no id
          set runloop False
        ]
        set rest-seedpool (rest-seedpool + (item id seedpool))
        set id id + 1
        if (id >= maxpft) [ set runloop False]
      ]

      ;;neue Art generieren
      let new-species (item new-species-no specieslist)
      sprout-trees 1 [set-species new-species]
    ]
  ]
end


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; disturbances ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
to disturbance
  let x random-pxcor
  let y random-pycor
  set-current-plot-pen "pen-1"
  set-plot-pen-mode 1
  ;;old: if (random 100 < disturbance-chance)

  let zufall random 100 / 100
  ifelse (zufall < disturbance-chance)
  [
     ask patch x y
     [
          plot 100
        	;;old: ask trees in-radius random-exponential disturbance-range
          let current random max-disturbance-radius
          let disturbance-range min-disturbance-radius + max(list 0 (current - min-disturbance-radius))
         ask trees in-radius disturbance-range
         [
          	  set pcolor magenta
          	  ;;old: if (random 100 < random-poisson disturbance-severity)[
          	  set status "dead"
          	  set no_died no_died + 1
          	  ;;old: ]
         ]
      ]
  ]
  [ ; else
      plot 0
  ]
  set-current-plot-pen "pen-0"
  set-plot-pen-mode 0
end


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; update tree status  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to update-status
  ask trees
  [
    if (status = "dead") [
      set former_species species
      set former_shadetolerance shadetolerance
      ask patch-here [ set pcolor 23 ]
      die
    ]
    if (status = "juvenile") [set status "alive"]
  ]
end


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Plotting Procedures ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to create-plot
  ;set-current-plot "Percentage monodom trees"
  clear-plot
end
@#$#@#$#@
GRAPHICS-WINDOW
424
10
932
519
-1
-1
5.0
1
10
1
1
1
0
1
1
1
0
99
0
99
1
1
1
ticks
30.0

BUTTON
38
62
103
95
NIL
go
T
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

BUTTON
38
27
103
60
NIL
setup
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

PLOT
981
14
1375
316
percentage Monodom
time
%
0.0
50.0
0.0
50.0
true
false
"" ""
PENS
"pen-0" 1.0 0 -16777216 true "" "plot count trees with [species = \"pft1\"] * 100 / (count trees)"
"pen-1" 1.0 0 -7500403 true "" ""

SLIDER
205
12
405
45
WorldSize
WorldSize
10
256
100.0
2
1
NIL
HORIZONTAL

SLIDER
203
54
375
87
seedmass
seedmass
8000
20000
10000.0
10
1
NIL
HORIZONTAL

SLIDER
37
129
232
162
disturbance-chance
disturbance-chance
0
0.2
0.1
0.001
1
NIL
HORIZONTAL

SLIDER
37
174
279
207
max-disturbance-radius
max-disturbance-radius
0
100
10.0
1
1
NIL
HORIZONTAL

SLIDER
37
216
260
249
min-disturbance-radius
min-disturbance-radius
0
20
10.0
1
1
NIL
HORIZONTAL

SWITCH
146
336
317
369
global-dispersal
global-dispersal
1
1
-1000

@#$#@#$#@
@#$#@#$#@
default
true
0
Polygon -7500403 true true 150 5 40 250 150 205 260 250

airplane
true
0
Polygon -7500403 true true 150 0 135 15 120 60 120 105 15 165 15 195 120 180 135 240 105 270 120 285 150 270 180 285 210 270 165 240 180 180 285 195 285 165 180 105 180 60 165 15

arrow
true
0
Polygon -7500403 true true 150 0 0 150 105 150 105 293 195 293 195 150 300 150

box
false
0
Polygon -7500403 true true 150 285 285 225 285 75 150 135
Polygon -7500403 true true 150 135 15 75 150 15 285 75
Polygon -7500403 true true 15 75 15 225 150 285 150 135
Line -16777216 false 150 285 150 135
Line -16777216 false 150 135 15 75
Line -16777216 false 150 135 285 75

bug
true
0
Circle -7500403 true true 96 182 108
Circle -7500403 true true 110 127 80
Circle -7500403 true true 110 75 80
Line -7500403 true 150 100 80 30
Line -7500403 true 150 100 220 30

butterfly
true
0
Polygon -7500403 true true 150 165 209 199 225 225 225 255 195 270 165 255 150 240
Polygon -7500403 true true 150 165 89 198 75 225 75 255 105 270 135 255 150 240
Polygon -7500403 true true 139 148 100 105 55 90 25 90 10 105 10 135 25 180 40 195 85 194 139 163
Polygon -7500403 true true 162 150 200 105 245 90 275 90 290 105 290 135 275 180 260 195 215 195 162 165
Polygon -16777216 true false 150 255 135 225 120 150 135 120 150 105 165 120 180 150 165 225
Circle -16777216 true false 135 90 30
Line -16777216 false 150 105 195 60
Line -16777216 false 150 105 105 60

car
false
0
Polygon -7500403 true true 300 180 279 164 261 144 240 135 226 132 213 106 203 84 185 63 159 50 135 50 75 60 0 150 0 165 0 225 300 225 300 180
Circle -16777216 true false 180 180 90
Circle -16777216 true false 30 180 90
Polygon -16777216 true false 162 80 132 78 134 135 209 135 194 105 189 96 180 89
Circle -7500403 true true 47 195 58
Circle -7500403 true true 195 195 58

circle
false
0
Circle -7500403 true true 0 0 300

circle 2
false
0
Circle -7500403 true true 0 0 300
Circle -16777216 true false 30 30 240

cow
false
0
Polygon -7500403 true true 200 193 197 249 179 249 177 196 166 187 140 189 93 191 78 179 72 211 49 209 48 181 37 149 25 120 25 89 45 72 103 84 179 75 198 76 252 64 272 81 293 103 285 121 255 121 242 118 224 167
Polygon -7500403 true true 73 210 86 251 62 249 48 208
Polygon -7500403 true true 25 114 16 195 9 204 23 213 25 200 39 123

cylinder
false
0
Circle -7500403 true true 0 0 300

dot
false
0
Circle -7500403 true true 90 90 120

face happy
false
0
Circle -7500403 true true 8 8 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Polygon -16777216 true false 150 255 90 239 62 213 47 191 67 179 90 203 109 218 150 225 192 218 210 203 227 181 251 194 236 217 212 240

face neutral
false
0
Circle -7500403 true true 8 7 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Rectangle -16777216 true false 60 195 240 225

face sad
false
0
Circle -7500403 true true 8 8 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Polygon -16777216 true false 150 168 90 184 62 210 47 232 67 244 90 220 109 205 150 198 192 205 210 220 227 242 251 229 236 206 212 183

fish
false
0
Polygon -1 true false 44 131 21 87 15 86 0 120 15 150 0 180 13 214 20 212 45 166
Polygon -1 true false 135 195 119 235 95 218 76 210 46 204 60 165
Polygon -1 true false 75 45 83 77 71 103 86 114 166 78 135 60
Polygon -7500403 true true 30 136 151 77 226 81 280 119 292 146 292 160 287 170 270 195 195 210 151 212 30 166
Circle -16777216 true false 215 106 30

flag
false
0
Rectangle -7500403 true true 60 15 75 300
Polygon -7500403 true true 90 150 270 90 90 30
Line -7500403 true 75 135 90 135
Line -7500403 true 75 45 90 45

flower
false
0
Polygon -10899396 true false 135 120 165 165 180 210 180 240 150 300 165 300 195 240 195 195 165 135
Circle -7500403 true true 85 132 38
Circle -7500403 true true 130 147 38
Circle -7500403 true true 192 85 38
Circle -7500403 true true 85 40 38
Circle -7500403 true true 177 40 38
Circle -7500403 true true 177 132 38
Circle -7500403 true true 70 85 38
Circle -7500403 true true 130 25 38
Circle -7500403 true true 96 51 108
Circle -16777216 true false 113 68 74
Polygon -10899396 true false 189 233 219 188 249 173 279 188 234 218
Polygon -10899396 true false 180 255 150 210 105 210 75 240 135 240

house
false
0
Rectangle -7500403 true true 45 120 255 285
Rectangle -16777216 true false 120 210 180 285
Polygon -7500403 true true 15 120 150 15 285 120
Line -16777216 false 30 120 270 120

leaf
false
0
Polygon -7500403 true true 150 210 135 195 120 210 60 210 30 195 60 180 60 165 15 135 30 120 15 105 40 104 45 90 60 90 90 105 105 120 120 120 105 60 120 60 135 30 150 15 165 30 180 60 195 60 180 120 195 120 210 105 240 90 255 90 263 104 285 105 270 120 285 135 240 165 240 180 270 195 240 210 180 210 165 195
Polygon -7500403 true true 135 195 135 240 120 255 105 255 105 285 135 285 165 240 165 195

line
true
0
Line -7500403 true 150 0 150 300

line half
true
0
Line -7500403 true 150 0 150 150

pentagon
false
0
Polygon -7500403 true true 150 15 15 120 60 285 240 285 285 120

person
false
0
Circle -7500403 true true 110 5 80
Polygon -7500403 true true 105 90 120 195 90 285 105 300 135 300 150 225 165 300 195 300 210 285 180 195 195 90
Rectangle -7500403 true true 127 79 172 94
Polygon -7500403 true true 195 90 240 150 225 180 165 105
Polygon -7500403 true true 105 90 60 150 75 180 135 105

plant
false
0
Rectangle -7500403 true true 135 90 165 300
Polygon -7500403 true true 135 255 90 210 45 195 75 255 135 285
Polygon -7500403 true true 165 255 210 210 255 195 225 255 165 285
Polygon -7500403 true true 135 180 90 135 45 120 75 180 135 210
Polygon -7500403 true true 165 180 165 210 225 180 255 120 210 135
Polygon -7500403 true true 135 105 90 60 45 45 75 105 135 135
Polygon -7500403 true true 165 105 165 135 225 105 255 45 210 60
Polygon -7500403 true true 135 90 120 45 150 15 180 45 165 90

sheep
false
15
Circle -1 true true 203 65 88
Circle -1 true true 70 65 162
Circle -1 true true 150 105 120
Polygon -7500403 true false 218 120 240 165 255 165 278 120
Circle -7500403 true false 214 72 67
Rectangle -1 true true 164 223 179 298
Polygon -1 true true 45 285 30 285 30 240 15 195 45 210
Circle -1 true true 3 83 150
Rectangle -1 true true 65 221 80 296
Polygon -1 true true 195 285 210 285 210 240 240 210 195 210
Polygon -7500403 true false 276 85 285 105 302 99 294 83
Polygon -7500403 true false 219 85 210 105 193 99 201 83

square
false
0
Rectangle -7500403 true true 30 30 270 270

square 2
false
0
Rectangle -7500403 true true 30 30 270 270
Rectangle -16777216 true false 60 60 240 240

star
false
0
Polygon -7500403 true true 151 1 185 108 298 108 207 175 242 282 151 216 59 282 94 175 3 108 116 108

target
false
0
Circle -7500403 true true 0 0 300
Circle -16777216 true false 30 30 240
Circle -7500403 true true 60 60 180
Circle -16777216 true false 90 90 120
Circle -7500403 true true 120 120 60

tree
false
0
Circle -7500403 true true 118 3 94
Rectangle -6459832 true false 120 195 180 300
Circle -7500403 true true 65 21 108
Circle -7500403 true true 116 41 127
Circle -7500403 true true 45 90 120
Circle -7500403 true true 104 74 152

triangle
false
0
Polygon -7500403 true true 150 30 15 255 285 255

triangle 2
false
0
Polygon -7500403 true true 150 30 15 255 285 255
Polygon -16777216 true false 151 99 225 223 75 224

truck
false
0
Rectangle -7500403 true true 4 45 195 187
Polygon -7500403 true true 296 193 296 150 259 134 244 104 208 104 207 194
Rectangle -1 true false 195 60 195 105
Polygon -16777216 true false 238 112 252 141 219 141 218 112
Circle -16777216 true false 234 174 42
Rectangle -7500403 true true 181 185 214 194
Circle -16777216 true false 144 174 42
Circle -16777216 true false 24 174 42
Circle -7500403 false true 24 174 42
Circle -7500403 false true 144 174 42
Circle -7500403 false true 234 174 42

turtle
true
0
Polygon -10899396 true false 215 204 240 233 246 254 228 266 215 252 193 210
Polygon -10899396 true false 195 90 225 75 245 75 260 89 269 108 261 124 240 105 225 105 210 105
Polygon -10899396 true false 105 90 75 75 55 75 40 89 31 108 39 124 60 105 75 105 90 105
Polygon -10899396 true false 132 85 134 64 107 51 108 17 150 2 192 18 192 52 169 65 172 87
Polygon -10899396 true false 85 204 60 233 54 254 72 266 85 252 107 210
Polygon -7500403 true true 119 75 179 75 209 101 224 135 220 225 175 261 128 261 81 224 74 135 88 99

wheel
false
0
Circle -7500403 true true 3 3 294
Circle -16777216 true false 30 30 240
Line -7500403 true 150 285 150 15
Line -7500403 true 15 150 285 150
Circle -7500403 true true 120 120 60
Line -7500403 true 216 40 79 269
Line -7500403 true 40 84 269 221
Line -7500403 true 40 216 269 79
Line -7500403 true 84 40 221 269

wolf
false
0
Polygon -16777216 true false 253 133 245 131 245 133
Polygon -7500403 true true 2 194 13 197 30 191 38 193 38 205 20 226 20 257 27 265 38 266 40 260 31 253 31 230 60 206 68 198 75 209 66 228 65 243 82 261 84 268 100 267 103 261 77 239 79 231 100 207 98 196 119 201 143 202 160 195 166 210 172 213 173 238 167 251 160 248 154 265 169 264 178 247 186 240 198 260 200 271 217 271 219 262 207 258 195 230 192 198 210 184 227 164 242 144 259 145 284 151 277 141 293 140 299 134 297 127 273 119 270 105
Polygon -7500403 true true -1 195 14 180 36 166 40 153 53 140 82 131 134 133 159 126 188 115 227 108 236 102 238 98 268 86 269 92 281 87 269 103 269 113

x
false
0
Polygon -7500403 true true 270 75 225 30 30 225 75 270
Polygon -7500403 true true 30 75 75 30 270 225 225 270
@#$#@#$#@
NetLogo 6.1.1
@#$#@#$#@
@#$#@#$#@
@#$#@#$#@
<experiments>
  <experiment name="disturbance_experiment" repetitions="200" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>go</go>
    <timeLimit steps="5000"/>
    <exitCondition>count turtles with [species = "pft1"] = 0</exitCondition>
    <metric>count turtles with [species = "pft1"] / count turtles</metric>
    <enumeratedValueSet variable="disturbance-chance">
      <value value="0"/>
      <value value="0.01"/>
      <value value="0.02"/>
      <value value="0.03"/>
      <value value="0.04"/>
      <value value="0.05"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-disturbance-radius">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-disturbance-radius">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="seedmass">
      <value value="98000"/>
    </enumeratedValueSet>
  </experiment>
</experiments>
@#$#@#$#@
@#$#@#$#@
default
0.0
-0.2 0 0.0 1.0
0.0 1 1.0 0.0
0.2 0 0.0 1.0
link direction
true
0
Line -7500403 true 150 150 90 180
Line -7500403 true 150 150 210 180
@#$#@#$#@
0
@#$#@#$#@
