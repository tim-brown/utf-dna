(ns utf8-dna.core
  "Implementation of [UTF-DNA: A Text Encoding for DNA Sequences](https://utf-dna.com/)"
  (:require [clojure.core.match :as m]))

(defn string->code-points
  "returns a sequence of the Unicode code points of s (optionaly starting at character i)"
  ([s i]
   (when (< i (count s))
     (let [cp (.codePointAt s i)]
       (lazy-seq (cons cp (string->code-points s (+ i (Character/charCount cp))))))))
  ([s] (string->code-points s 0)))

(def A "nucleiotide mapping for A" 0)
(def C "nucleiotide mapping for C" 1)
(def G "nucleiotide mapping for G" 2)
(def T "nucleiotide mapping for T" 3)

(def nucleotide-char-mapping {\A A \C C \G G \T T})

(defn utf8-dna-encode-code-point
  "encodes UNICODE code-point cp into a UTF8-DNA sequence of quaternary digits"
  [cp]
  (letfn [(nth-quaternary [n q] (bit-and 0x3 (bit-shift-right n (* q 2))))
          (T-codons [x & low-qs] (mapcat #(vector T (nth-quaternary x (inc %)) (nth-quaternary x %)) low-qs))]
    (cond
      (<= cp 0xf) (list A (nth-quaternary cp 1) (nth-quaternary cp 0))
      (<= cp 0xff) (list* C (nth-quaternary cp 3) (nth-quaternary cp 2) (T-codons cp 0))
      (<= cp 0x3ff) (list* G A (nth-quaternary cp 4) (T-codons cp 2 0))
      (<= cp 0x3fff) (list* G C (nth-quaternary cp 6) (T-codons cp 4 2 0))
      (<= cp 0x3ffff) (list* G G (nth-quaternary cp 8) (T-codons cp 6 4 2 0))
      (<= cp 0x3fffff) (list* G T (nth-quaternary cp 10) (T-codons cp 8 6 4 2 0)))))

(defn string->utf-dna
  "returns a lazy sequence of ACGT characters encoding s in UTF8-DNA"
  [s]
  (->> (string->code-points s)
      (mapcat utf8-dna-encode-code-point)
      (map {0 \A 1 \C 2 \G 3 \T})))

(defn dna->code-points
  "returns a lazy sequence of UNICODE code points represented by the nucleioides sequence"
  [nucleotides]
  (letfn [(reconstitute [& bits] (reduce #(+ (* %1 4) %2) 0 bits))]
    (when (seq nucleotides)
      (m/match [nucleotides]
        [([0 Z z & r] :seq)] (lazy-seq (cons (reconstitute Z z) (dna->code-points r)))
        [([1 Y y, 3 Z z & r] :seq)] (lazy-seq (cons (reconstitute Y y Z z) (dna->code-points r)))
        [([2 0 x, 3 Y y, 3 Z z & r] :seq)] (lazy-seq (cons (reconstitute x Y y Z z) (dna->code-points r)))
        [([2 1 w, 3 X x, 3 Y y, 3 Z z & r] :seq)] (lazy-seq (cons (reconstitute w X x Y y Z z) (dna->code-points r)))
        [([2 2 v, 3 W w, 3 X x, 3 Y y, 3 Z z & r] :seq)] (lazy-seq (cons (reconstitute v W w X x Y y Z z) (dna->code-points r)))
        [([2 3 u, 3 V v, 3 W w, 3 X x, 3 Y y, 3 Z z & r] :seq)] (lazy-seq (cons (reconstitute u V v W w X x Y y Z z) (dna->code-points r)))))))

(defn utf8-dna->string
  "returns a lazy sequence of characters encoded by the dna string"
  [dna]
  (->> (keep nucleotide-char-mapping dna) dna->code-points (mapcat Character/toChars)))

(comment
  (utf8-dna->string "ACG") ; single A..
  (utf8-dna->string "CTATAA") ; single C..T.. (probably your most common case)
  (utf8-dna->string "CTATAACAGTAA") ; sequential C..T..{2}
  (utf8-dna->string "CTATAACAGTAACCTTCG") ; longer C...T..{>2}
  (utf8-dna->string "GGC TTT TCG TAA TGG") ; GG.T..T..T..T..
  (utf8-dna->string "
CTA TAA CAG TAA CCT TCG CCG TTT CCT TCA CCT
TAG CCG TCC CAG TAA CCT TAT CCG TAC CCG TTG
CCT TCA CTG TGC CAG TAA GGC TTT TCG TAA TGG
  ")
  (utf8-dna->string "
CTATAACAGTAACCTTCGCCGTTTCCTTCACCT
TAGCCGTCCCAGTAACCTTATCCGTACCCGTTG
CCTTCACTGTGCCAGTAAGGCTTTTCGTAATGG
                    ")
  (string->code-points "À votre santé 😊")
  (utf8-dna->string (string->utf-dna "À votre santé 😊"))
  (string->utf-dna "À votre santé 😊"))