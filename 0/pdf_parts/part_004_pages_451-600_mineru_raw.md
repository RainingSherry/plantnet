It might therefore seem best to ignore the other, rather subtle representations and just look at these ones. But it turns out to be easier to study all representations and only later ask which ones occur inside $L ^ { 2 } ( X )$ . For $\mathrm { S L } _ { 2 } ( \mathbb { R } )$ , the representations we have just constructed (which were subquotients of $W _ { s } ^ { \pm } )$ exhaust all the irreducible representations,7 and there is a Plancherel formula for $L ^ { 2 } ( \mathbb { R } ^ { 2 } )$ that tells us which ones appear in $L ^ { 2 } ( \mathbb { R } ^ { 2 } )$ and with what multiplicity:

$$
L ^ {2} (\mathbb {R} ^ {2}) = \int_ {- \infty} ^ {\infty} W _ {- 1 + \mathrm{i} t} \mathrm{e} ^ {\mathrm{i} t} \mathrm{d} t.
$$

To summarize: if G is not compact, then we can no longer take averages over G. This has various consequences:

Representations occur in continuous families. The decomposition of $L ^ { 2 } ( X )$ takes the form of a direct integral, not a direct sum.

Representations do not split up into a direct sum of irreducibles. Even when a representation admits a finite composition series, as with the action of $\mathrm { S L } ^ { 2 } ( \mathbb { R } )$ on $W _ { s } ^ { \pm }$ , it need not split up into a direct sum. So to describe all representations we need to do more than just describe the irreducibles—we also need to describe the glue that holds them together.

So far, the theory of representations of a noncompact group G seems to have none of the pleasant features of the compact case. But one thing does survive: there is still an analogue of the theorem that the character table is square. Indeed, we can still define characters in terms of the traces of group elements. But now we must be careful, since the irreducible representation may be on an infinite-dimensional vector space, so that its trace cannot be defined so easily. In fact, characters are not functions on G, but only distributions [III.18]. The character of a representation determines the semisimplification of a representation ρ: that is, it tells us which irreducible representations are part of $\rho ,$ but not how they are glued together.8

These phenomena were discovered by Harish-Chandra in the 1950s in an extraordinary series of works that completely described the representation theory of Lie groups such as the ones we have discussed (the precise condition is that they should be real and reductive— a concept that will be explained later in this article) and the generalizations of classical theorems of Fourier analysis to this setting.9

Independently and slightly earlier, Brauer had investigated the representation theory of finite groups on finite-dimensional vector spaces over fields of characteristic p. Here, too, reducible representations need not decompose as direct sums, though in this case the problem is not lack of compactness (obviously, since everything is finite) but an inability to average over the group: we would like to divide by |G|, but often this is zero. A simple example that illustrates this is the action of $\mathbb { Z } / p \mathbb { Z }$ on the space $\mathbb { F } _ { p } ^ { 2 }$ that takes x to the $2 \times 2$ matrix ${ \bigl ( } { } _ { 1 } ^ { 1 } \mathbf { \Sigma } _ { 0 } ^ { x } { \bigr ) }$ . This is reducible, since the column vector $\bigl ( \begin{array} { l } { 1 } \\ { 0 } \end{array} \bigr )$ is fixed by the action, and therefore generates an invariant subspace. However, if one could decompose the action, then the matrices $\big (  _ { 1 0 } ^ { 1 x } \big )$ would all be diagonalizable, which they are not.

It is possible for there to be infinitely many indecomposable representations, which again may vary in families. However, as before, there are only finitely many irreducible representations, so there is some chance of a “character table is square” theorem in which the rows of the square are parametrized by characters of irreducible representations. Brauer proved just such a theorem, pairing the characters with p-semisimple conjugacy classes in G: that is, conjugacy classes of elements whose order is not divisible by p.

We will draw two crude morals from the work of Harish-Chandra and of Brauer. The first is that the category of representations of a group is always a reasonable object, but when the representations are infinite dimensional it requires serious technical work to set it up. Objects in this category do not necessarily decompose as a direct sum of irreducibles (one says that the category is not semisimple), and can occur in infinite families, but irreducible objects pair off in some precise way with certain “diagonalizable” conjugacy classes in the group—there is always some kind of analogue of “the character table is square” theorem.

It turns out that when we consider representations in more general contexts—Lie algebras acting on vector spaces, quantum groups, p-adic groups on infinitedimensional complex or p-adic vector spaces, etc.— these qualitative features stay the same.

The second moral is that we should always hope for some “non-Abelian Fourier transform”: that is, a set that parametrizes irreducible representations and a description of the character values in terms of this set.

In the case of real reductive groups Harish-Chandra’s work provides such an answer, generalizing the Weyl character formula for compact groups; for arbitrary groups no such answer is known. For special classes of groups, there are partially successful general principles (the orbit method, Broué’s conjecture), of which the deepest are the extraordinary circle of conjectures known as the Langlands program, which we shall discuss later.

# 5 Interlude: The Philosophical Lessons of “The Character Table Is Square”

Our basic theorem (“the character table is square”) tells us to expect that the category of all irreducible representations of G is interesting when the conjugacyclass structure of G is in some way under control. We will finish this essay by explaining a remarkable family of examples of such groups—the rational points of reductive algebraic groups—and their conjectured representation theory, which is described by the Langlands program.

An affine algebraic group is a subgroup of some group ${ \mathrm { G L } } _ { n }$ that is defined by polynomial equations in the matrix coefficients. For example, the determinant of a matrix is a polynomial in the matrix coefficients, so the group $\mathrm { S L } _ { n } ,$ , which consists of all matrices in ${ \mathrm { G L } } _ { n }$ with determinant 1, is such a group. Another is ${ \mathrm { S O } } _ { n }$ , which is the set of matrices with determinant 1 that satisfy the equation $A A ^ { \mathrm { T } } = I .$ .

The above notation did not specify what sort of coefficients we were allowing for the matrices. That vagueness was deliberate. Given an algebraic group G and a field k, let us write G(k) for the group where the coefficients are taken to have values in k. For example, $\mathrm { S L } _ { n } ( \mathbb { F } _ { q } )$ is the set of n  n matrices with coefficients in the finite field $\mathbb { F } _ { q }$ and determinant 1. This group is finite, as is ${ \mathrm { S O } } _ { n } ( \mathbb { F } _ { q } )$ , while $\mathrm { S L } _ { n } ( \mathbb { R } )$ and $\mathrm { S O } _ { n } ( \mathbb { R } )$ are Lie groups. Moreover, $\mathrm { S O } _ { n } ( \mathbb { R } )$ is compact, while $\mathrm { S L } _ { n } ( \mathbb { R } )$ is not. So among affine algebraic groups over fields one already finds all three types of groups we have discussed: finite groups, compact Lie groups, and noncompact Lie groups.

We can think of $\mathrm { S L } _ { n } ( \mathbb { R } )$ as the set of matrices in $\mathrm { S L } _ { n } ( \mathbb { C } )$ that are equal to their complex conjugates. There is another involution on $\mathrm { S L } _ { n } ( \mathbb { C } )$ that is a sort of “twisted” form of complex conjugation, where we send a matrix A to the complex conjugate of $( A ^ { - 1 } ) ^ { \mathrm { T } }$ . The fixed points of this new involution (that ${ \mathrm { i } } \mathbf { s } ,$ the determinant-1 matrices A such that A equals the complex conjugate of $( A ^ { - 1 } ) ^ { \mathrm { T } } )$ form a group called $\mathrm { S U } _ { n } ( \mathbb { R } )$ . This is also called a real form of $\mathrm { S L } _ { n } ( \mathbb { C } ) , ^ { 1 0 }$ and it is compact.

The groups $\mathrm { S L } _ { n } ( \mathbb { F } _ { q } )$ and ${ \mathrm { S O } } _ { n } ( \mathbb { F } _ { q } )$ are almost simple groups;11 the classification of finite simple groups tells us, mysteriously, that all but twenty-six of the finite simple groups are of this form. A much, much easier theorem tells us that the connected compact groups are also of this form.

Now, given an algebraic group $G ,$ we can also consider the instances $G ( \mathbb { Q } _ { p } )$ , where $\mathbb { Q } _ { p }$ is the field of p-adic numbers, and also G(Q). For that matter, we may consider G(k) for any other field $k ,$ such as the function field of an algebraic variety [V.30]. The lesson of section 4 is that we may hope for all of these many groups to have a good representation theory, but that to obtain it there will be serious “analytic” or “arithmetic” difficulties to overcome, which will depend strongly on the properties of the field k.

Lest the reader adopt too optimistic a viewpoint, we point out that not every affine algebraic group has a nice conjugacy-class structure. For example, let $V _ { n }$ be the set of upper triangular matrices in ${ \mathrm { G L } } _ { n }$ with 1s along the diagonal, and let k be $\mathbb { F } _ { q } .$ . For large n, the conjugacy classes in $V _ { n } ( \mathbb { F } _ { q } )$ form large and complex families: to parametrize them sensibly one needs more than n parameters (in other words, they belong to families of dimension greater than n, in an appropriate sense), and it is not in fact known how to parametrize them even for a smallish value of $n ,$ such as 11. (It is not obvious that this is a “good” question though.)

More generally, solvable groups tend to have horrible conjugacy-class structure, even when the groups themselves are “sensible.” So we might expect their representation theory to be similarly horrible. The best we can hope for is a result that describes the entries of the character table in terms of this horrible structure— some kind of non-Abelian Fourier integral. For certain p-groups Kirillov found such a result in the 1960s, as an example of the “orbit method,” but the general result is not yet known.

On the other hand, groups that are similar to connected compact groups do have a nice conjugacy-class structure: in particular, finite simple groups do. An algebraic group is called reductive if G(C) has a compact real form. So, for instance, $\mathrm { S L } _ { n }$ is reductive by the existence of the real form $\mathrm { S U } _ { n } ( \mathbb { R } )$ . The groups ${ \mathrm { G L } } _ { n }$ and $\mathrm { S O } _ { n }$ are also reductive, but $V _ { n }$ is not.12

Let us examine the conjugacy classes in the group $\operatorname { S U } _ { n }$ . Every matrix in $\mathrm { S U } _ { n } ( \mathbb { R } )$ can be diagonalized, and two conjugate matrices have the same eigenvalues, up to reordering. Conversely, any two matrices in $\mathrm { S U } _ { n } ( \mathbb { R } )$ with the same eigenvalues are conjugate. Therefore, the conjugacy classes are parametrized by the quotient of the subgroup of all diagonal matrices by the action of $S _ { n }$ that permutes the entries.

This example can be generalized. Any compact connected group has a maximal torus T , that is, a maximal subgroup isomorphic to a product of circles. (In the previous example it was the subgroup of diagonal matrices.) Any two maximal tori are conjugate in $G ,$ and any conjugacy class in G intersects T in a unique W -orbit on T , where W is the Weyl group, the finite group $N ( T ) / T$ (where N(T ) is the normalizer of T ).

The description of conjugacy classes in $G ( { \bar { k } } )$ , for an algebraically closed field k¯, is only a little more complicated. Any element $g \in G ( \bar { k } )$ admits a jordan decomposition [III.43]: it can be written as $g \ = \ s u \ = \ u s ,$ , where s is conjugate to an element of $T ( \bar { k } )$ and u is unipotent when considered as an element of ${ \mathrm { G L } } _ { n } ( { \bar { k } } )$ . (A matrix A is unipotent if some power of $A - I$ is zero.) Unipotent elements never intersect compact subgroups. When $G = { \mathrm { G L } } _ { n }$ this is the usual Jordan decomposition; conjugacy classes of unipotent elements are parametrized by partitions of n, which, as we mentioned in section 2, are precisely the conjugacy classes of $W = S _ { n } ,$ . For general reductive groups, unipotent conjugacy classes are again almost the same thing as conjugacy classes in $W . ^ { 1 3 }$ In particular, there are finitely many, independent of k¯.

Finally, when k is not algebraically closed, one describes conjugacy classes by a kind of Galois descent; for example, in ${ \mathrm { G L } } _ { n } ( k )$ , semisimple classes are still determined by their characteristic polynomial, but the fact that this polynomial has coefficients in k constrains the possible conjugacy classes.

The point of describing the conjugacy-class structure in such detail is to describe the representation theory in analogous terms. A crude feature of the conjugacyclass structure is the way it decouples the field k from finite combinatorial data that is attached to G but independent of k—things like W , the lattice defining $T ,$ roots, and weights.

The “philosophy” suggested by the theorem that the character table is square suggests that the representation theory should also admit such a decoupling: it should be built out of the representation theory of $k ^ { * }$ , which is the analogue of the circle, and out of the combinatorial structure of $G ( { \bar { k } } )$ (such as the finite groups W ). Moreover, representations should have a “Jordan decomposition”:14 the “unipotent” representations should have some kind of combinatorial complexity but little dependence on k, and compact groups should have no unipotent representations.

The Langlands program provides a description along the lines laid out above, but it goes beyond any of the results we have suggested in that it also describes the entries of the character table. Thus, for this class of examples, it gives us (conjecturally) the hoped-for “non-Abelian Fourier transform.”

# 6 Coda: The Langlands Program

And so we conclude by just hinting at statements. If G(k) is a reductive group, we want to describe an appropriate category of representations for G(k), or at least the character table, which we may think of as a “semisimplification” of that category.

Even when k is finite, it is too much to hope that conjugacy classes in G(k) parametrize irreducible representations. But something not so far off is conjectured, as follows.

To a reductive group G over an algebraically closed field, Langlands attaches another reductive group $^ \mathrm { L } G ,$ the Langlands dual, and conjectures that representations of G(k) will be parametrized by conjugacy classes in ${ } ^ { \mathrm { L } } G ( \mathbb { C } ) . { } ^ { 1 5 }$ However, these are not conjugacy classes of elements of $\operatorname { L } G ( \mathbb { C } )$ , as before, but of homomorphisms from the Galois group of $k ~ { \mathrm { t o } } ~ ^ { \mathrm { L } } G$ . The Langlands dual was originally defined in a combinatorial manner, but there is now a conceptual definition. A few examples of pairs $( G , ^ { \mathrm { L } } G )$ are $( { \mathrm { G L } } _ { n } , { \mathrm { G L } } _ { n } )$ , $( \mathrm { S O } _ { 2 n + 1 } , \mathrm { S p } _ { 2 n } )$ , and $( \mathrm { S L } _ { n } , \mathrm { P G L } _ { n } )$ .

In this way the Langlands program describes the representation theory as built out of the structure of G and the arithmetic of k.

Although this description indicates the flavor of the conjectures, it is not quite correct as stated. For instance, one has to modify the Galois group16 in such a way that the correspondence is true for the group $\operatorname { G L } _ { 1 } ( k ) = k ^ { * }$ . When $k = \mathbb { R } ,$ , we get the representation theory of R∗ (or its compact form $S ^ { 1 } )$ , which is Fourier analysis; on the other hand, when k is a p-adic local field, the representation theory of k∗ is described by local class field theory. We already see an extraordinary aspect of the Langlands program: it precisely unifies and generalizes harmonic analysis and number theory.

The most compelling versions of the Langlands program are “equivalences of derived categories” between the category of representations and certain geometric objects on the spaces of Langlands parameters. These conjectural statements are the hoped-for Fourier transforms.

Though much progress has been made, a large part of the Langlands program remains to be proved. For finite reductive groups, slightly weaker statements have been proved, mostly by Lusztig. As all but twenty-six of the finite simple groups arise from reductive groups, and as the sporadic groups have had their character tables computed individually, this work already determines the character tables of all the finite simple groups.

For groups over R, the work of Harish-Chandra and later authors again confirms the conjectures. But for other fields, only fragmentary theorems have been proved. There is much still to be done.

# Further Reading

A nice introductory text on representation theory is Alperin’s Local Representation Theory (Cambridge University Press, Cambridge, 1993). As for the Langlands program, the 1979 American Mathematical Society volume titled Automorphic Forms, Representations, and L-functions (but universally known as “The Corvallis Proceedings”) is more advanced, and as good a place to start as any.

# IV.10 Geometric and Combinatorial Group Theory

Martin R. Bridson

# 1 What Are Combinatorial and Geometric Group Theory?

Groups and geometry are ubiquitous in mathematics, groups because the symmetries (or automorphisms [I.3 §4.1]) of any mathematical object in any context form a group and geometry because it allows one to think intuitively about abstract problems and to organize families of objects into spaces from which one may gain some global insight.

The purpose of this article is to introduce the reader to the study of infinite, discrete groups. I shall discuss both the combinatorial approach to the subject that held sway for much of the twentieth century and the more geometric perspective that has led to an enormous flowering of the subject in the last twenty years. I hope to convince the reader that the study of groups is a concern for all of mathematics rather than something that belongs particularly to the domain of algebra.

The principal focus of geometric group theory is the interaction of geometry/topology and group theory, through group actions and through suitable translations of geometric concepts into group theory. One wants to develop and exploit this interaction for the benefit of both geometry/topology and group theory. And, in keeping with our assertion that groups are important throughout mathematics, one hopes to illuminate and solve problems from elsewhere in mathematics by encoding them as problems in group theory.

Geometric group theory acquired a distinct identity in the late 1980s but many of its principal ideas have their roots in the end of the nineteenth century. At that time, low-dimensional topology and combinatorial group theory emerged entwined. Roughly speaking, combinatorial group theory is the study of groups defined in terms of presentations, that is, by means of generators and relations. In order to follow the rest of this introduction the reader must first understand what these terms mean. Since their definitions would require an unacceptably long break in the flow of our discussion, I will postpone them to the next section, but I strongly advise the reader who is unfamiliar with the meaning of the expression $\boldsymbol { \Gamma } = \langle a _ { 1 } , \dots , a _ { n } \mid r _ { 1 } , \dots , r _ { m } \rangle$ to pause and read that section before continuing with this one.

The rough definition of combinatorial group theory just given misses the point that, like many parts of mathematics, it is a subject defined more by its core problems and its origins than by its fundamental definitions. The initial impetus for the subject came from the description of discrete groups of hyperbolic isometries and, most particularly, the discovery of the fundamental group [IV.6 §2] of a manifold [I.3 §6.9] by poincaré [VI.61] in 1895. The group-theoretic issues that emerged were brought into sharp focus by the work of Tietze and Dehn in the first decade of the twentieth century and drove much of combinatorial group theory for the remainder of the century.

Not all of the epoch-defining problems came from topology: other areas of mathematics threw up fundamental questions as well. Here are some of the forms they took: Does there exist a group of the following type? Which groups have the following property? What are the subgroups of …? Is the following group infinite? When can one determine the structure of a group from its finite quotients? In the sections that follow I shall attempt to illustrate the mathematical culture associated with questions of this kind, but let me immediately mention some easily stated but difficult classical problems. (i) Let G be a group that is finitely generated and suppose that there is some positive integer n such that $x ^ { n } = 1$ for every x in G. Must G be finite? (ii) Is there a finitely presented group Γ and a surjective homomorphism φ : Γ → Γ such that φ(γ) = 1 for some $\gamma \neq 1 2$ (iii) Does there exist a finitely presented, infinite, simple group [I.3 §3.3]? (iv) Is every countable group isomorphic to a subgroup of a finitely generated group, or even a finitely presented group?

The first of these questions was asked by Burnside in 1902 and the second by Hopf in connection with his study of degree-1 maps between manifolds. I shall present the answers to all four questions (in section 5) to illustrate an important aspect of both combinatorial and geometric group theory: one develops techniques that allow the construction of explicit groups with prescribed properties. Such constructions are of particular interest when they illustrate the diversity of possible phenomena in other branches of mathematics.

Another kind of question that raises basic issues in combinatorial group theory takes the form: Does there exist an algorithm to determine whether or not a group (or given elements of a group) has such-andsuch a property? For example, does there exist an algorithm that can take any finite presentation and decide in a finite number of steps whether or not the group presented is trivial? Questions of this type led to a profound and mutually beneficial interaction between group theory and logic, given full voice by the Higman embedding theorem, which we shall discuss in section 6. Moreover, via the conduit of combinatorial group theory, logic has influenced topology as well: one uses group-theoretic constructions to show, for example, that there is no algorithm to determine which pairs of compact triangulated manifolds are homeomorphic in dimensions 4 and above. This shows that certain kinds of classification results that have been obtained in two and three dimensions do not have higher-dimensional analogues.

One might reasonably regard combinatorial group theory as the attempt to develop algebraic techniques to solve the types of questions described above, and in the course of doing so to identify classes of groups that are worthy of particular study. This last point, the question of which groups deserve our attention, is tackled head-on in the final section of this article.

Some of the triumphs of combinatorial group theory are intrinsically combinatorial in nature, but many more have had their true nature revealed by the introduction of geometric techniques in the past twenty years. A fine example of this is the way in which Gromov’s insights have connected algorithmic problems in group theory to so-called filling problems in Riemannian geometry. Moreover, the power of geometric group theory is by no means confined to improving the techniques of combinatorial group theory: it naturally leads one to think about many other issues of fundamental importance. For example, it provides a context in which one can illuminate and vastly extend classical rigidity theorems [V.23], such as that of Mostow. The key to applications such as this is the idea that finitely generated groups can usefully be regarded as geometric objects in their own right. This idea has its origins in the work of cayley [VI.46] (1878) and Dehn (1905) but its full force was recognized and promoted by Gromov, starting in the 1980s. It is the key idea that underpins the later sections of this article.

# 2 Presenting Groups

How should one describe a group? An example will illustrate the standard way of doing so and give some idea of why it is often appropriate.

Consider the familiar tiling of the Euclidean plane by equilateral triangles. How might you describe the full group ${ \varGamma } _ { \Delta }$ of symmetries of this tiling, i.e., the rigid motions of the plane that send tiles to tiles? Let us focus on a single tile $T$ and a particular edge e of $T ,$ and use this to pick out three symmetries. The first, which we shall call $\alpha ,$ is the reflection of the plane in the line that contains e and the other two, β and $\gamma ,$ are the reflections in the lines that join the endpoints of e to the midpoints of the opposite edges in $T .$ With some effort one can convince oneself that every symmetry of the tiling can be obtained by performing these three operations repeatedly in a suitable order. One expresses this by saying that the set $\{ \alpha , \beta , \gamma \}$ generates the group $\Gamma _ { \Delta } .$

A further useful observation is that if one performs the operation α twice, the tiling is returned to its original position: that $\mathbf { i s } , \alpha ^ { 2 } = 1$ . Likewise, $\beta ^ { 2 } = \gamma ^ { 2 } = 1$ . One can also verify that $( \alpha \beta ) ^ { 6 } = ( \alpha \gamma ) ^ { 6 } = ( \beta \gamma ) ^ { 3 } = 1$ .

It turns out that the group $T _ { \Delta }$ is completely determined by these facts alone, a statement that we summarize by the notation

$$
\Gamma_ {\Delta} = \langle \alpha , \beta , \gamma | \alpha^ {2}, \beta^ {2}, \gamma^ {2}, (\alpha \beta) ^ {6}, (\alpha \gamma) ^ {6}, (\beta \gamma) ^ {3} \rangle .
$$

The aim of the rest of this section is to say in more detail what this means.

To begin with, notice that from the facts we are given we can deduce others: for example, bearing in mind that $\beta ^ { 2 } = \gamma ^ { 2 } = ( \beta \gamma ) ^ { 3 } = 1$ , we can show that

$$
(\gamma \beta) ^ {3} = (\gamma \beta) ^ {3} (\beta \gamma) ^ {3} = 1
$$

as well (where the last equality follows after repeatedly canceling pairs of the form ββ or γγ). We wish to convey the idea that in $T _ { \Delta }$ there are no relationships between the generators except those that follow from the facts above by this kind of argument.

Now let us try to say this more formally. We define a set of generators for a group Γ to be a subset $S \subset T$ such that every element of Γ is equal to some product of elements of S and their inverses. That is, every element can be written in the form $s _ { 1 } ^ { \varepsilon _ { 1 } } s _ { 2 } ^ { \varepsilon _ { 2 } } \ldots s _ { n } ^ { \varepsilon _ { n } }$ , where each $s _ { i }$ is an element of S and each $\varepsilon _ { i }$ is 1 $\mathrm { o r } - 1$ . We then call a product of this kind a relation if it is equal to the identity in Γ .

There is an awkward ambiguity here. When we talk about “the product” of some elements of Γ , it sounds as though we are referring to another element of Γ , but we certainly did not mean this at the end of the last paragraph: a relation is not the identity element of Γ but rather a string of symbols such as $a b ^ { - 1 } a ^ { - 1 } b c$ that yields the identity in Γ when you interpret $a , b ,$ , and c as generators in the set S. In order to be clear about this, it is useful to define another group, known as the free group F(S).

For concreteness we shall describe the free group with three generators, taking our set S to be $\{ a , b , c \}$ . A typical element is a “word” in the elements of S and their inverses, such as the expression ab $^ { - 1 } a ^ { - 1 } b c$ considered in the previous paragraph. However, we sometimes regard two words as the same: for instance, $a b c c ^ { - 1 } a c$ and $a b a b ^ { - 1 } b c$ are the same because they become identical when we cancel out the inverse pairs $c c ^ { - 1 }$ and $b ^ { - 1 } b$ . More formally, we define two such words to be equivalent and say that the elements of the free group are the equivalence classes [I.2 §2.3]. To multiply words together, we just concatenate them: for instance, the product of $a b ^ { - 1 }$ and bcca is $a b ^ { - 1 }$ bcca, which we can shorten to acca. The identity is the “empty word.” This is the free group on three generators a, b, and c. It should be clear how to generalize it to an arbitrary set S, though we shall continue to discuss the set $S = \{ a , b , c \}$ .

A more abstract way of characterizing the free group on a, $^ { b , }$ and c is to say that it has the following universal property: if G is any group and φ is any function from $S ~ = ~ \{ a , b , c \}$ to G, then there is a unique homomorphism Φ from F(S) to G that takes a to $\phi ( a )$ , b to φ(b), and c to φ(c). Indeed, if we want Φ to have these properties, then our definition is forced upon us: for example, $\phi ( a b ^ { - 1 } c a )$ will have to be $\phi ( a ) \phi ( b ) ^ { - 1 } \phi ( c ) \phi ( a )$ , by the definition of a homomorphism. So the uniqueness is obvious. The rough reason that this definition really does give rise to a well-defined homomorphism is that the only equations that are true in F(S) are ones that are true in all groups: in order for Φ not to be a homomorphism, one would need a relation to hold in F(S) that did not hold in $G ,$ but this is impossible.

Now let us return to our example ${ \cal { T } } _ { \Delta }$ . We would like to prove that it is (isomorphic to) the “freest” group with generators $\alpha , \beta ,$ and γ that satisfies the relations $\alpha ^ { 2 } ~ = ~ \beta ^ { 2 } ~ = ~ \gamma ^ { 2 } ~ = ~ ( \alpha \beta ) ^ { 6 } ~ = ~ ( \alpha \gamma ) ^ { 6 } ~ = ~ ( \beta \gamma ) ^ { 3 } ~ = ~ 1$ . But what exactly is this “freest” group that we are claiming is isomorphic to $I _ { \Delta } ?$

To avoid confusion about the meaning of $\alpha , \beta ,$ , and γ (are they elements of ${ \varGamma } _ { \Delta }$ or of the group that we are trying to construct that will turn out to be isomorphic to $I _ { \Delta } ? )$ we shall use the letters a, b, and c when we answer this question. Thus, we are trying to build the “freest” group with generators $a , b ,$ , and c that satisfies the relations $a ^ { 2 } ~ = ~ b ^ { 2 } ~ = ~ c ^ { 2 } ~ = ~ ( a b ) ^ { 6 } ~ = ~$ $( a c ) ^ { 6 } = ( b c ) ^ { 3 } = 1$ , which we denote by $G = \langle a , b , c \mid$ $a ^ { 2 } , b ^ { 2 } , c ^ { 2 } , ( a b ) ^ { 6 } , ( a c ) ^ { 6 } , ( b c ) ^ { 3 } \rangle$ .

There are two ways of going about this task. One is to imitate the above discussion of the free group itself, except that now we say that two words are equivalent if you can get from one to the other by inserting or deleting not just inverse pairs but also one of the words $a ^ { 2 } , b ^ { 2 } , c ^ { 2 } , ( a b ) ^ { 6 } , ( a c ) ^ { 6 } , \mathrm { o r } \ ( b c ) ^ { 3 }$ . For example, $a b ^ { 2 } c$ is equivalent to ac in this group. G is then defined to be the set of equivalence classes of words with the product coming from concatenation.

A neater way to obtain G is more conceptual and exploits the universal property of the free group. As G is to be generated by a, b, and $^ { c , }$ the universal property of the free group F(S) tells us that there will have to be a unique homomorphism Φ from $F ( S )$ to G such that $\phi ( a ) = a , \phi ( b ) = b ,$ and $\phi ( c ) = c ,$ . Moreover, we require that all of $a ^ { 2 } , b ^ { 2 } , c ^ { 2 } , ( a b ) ^ { 6 } , ( a c ) ^ { 6 }$ , and $( b c ) ^ { 3 }$ must map to the identity element in G. It follows that the kernel [I.3 §4.1] of Φ is a normal subgroup [I.3 §3.3] of F(S) that contains the set $R \ =$ $\{ a ^ { 2 } , b ^ { 2 } , c ^ { 2 } , ( a b ) ^ { 6 } , ( a c ) ^ { 6 } , ( b c ) ^ { 3 } \}$ . Let us write R for the smallest normal subgroup of $F ( S )$ that contains R (or equivalently the intersection of all normal subgroups of $F ( S )$ that contain R). Then there is a surjective homomorphism from the quotient [I.3 §3.3] $F ( S ) / \left. \left. R \right. \right.$ to any group that is generated by $a , b ,$ and c and satisfies the relations $a ^ { 2 } = b ^ { 2 } = c ^ { 2 } = ( a b ) ^ { 6 } =$ $( a c ) ^ { 6 } = ( b c ) ^ { 3 } = 1$ . This quotient itself is the group we are looking for: it is the largest group generated by $^ { a , }$ b, and c that satisfies the relations in R.

Our assertion about ${ \cal { T } } _ { \Delta }$ is that it is isomorphic to the group $G = \langle a , b , c \mid a ^ { 2 } , b ^ { 2 } , c ^ { 2 } , ( a b ) ^ { 6 } , ( a c ) ^ { 6 } , ( b c ) ^ { 3 } \rangle$ that we have just described (in two ways). More precisely, the map from $F ( S ) / \left. \left. R \right. \right.$ to ${ \varGamma } _ { \Delta }$ that takes a to $\alpha , b$ to $\beta ,$ , and c to γ is an isomorphism.

The above construction is very general. If we are given a group Γ , then a presentation of Γ is a set S that generates Γ , together with a set $R \subset F ( S )$ of relations, such that Γ is isomorphic to the quotient $F ( S ) / \left. \left. R \right. \right.$ . If both S and R are finite sets, one says that the presentation is finite. A group is finitely presented if it has a finite presentation.

We can also define presentations in the abstract, without mentioning a group Γ in advance: given any set S and any subset $R \subset F ( S )$ , we just define $\langle S \mid R \rangle$ to be the group $F ( S ) / \left. \left. R \right. \right.$ . This is the “freest” group generated by S that satisfies the relations in R: the only relations that hold in $\langle S \mid R \rangle$ are the ones that can be deduced from the relations R.

A psychological advantage of switching to this more abstract setting is that, whereas previously we began with a group Γ and asked how we might present it, we can now write down group presentations at will, starting with any set S and prescribing a set of words R in the symbols $S ^ { \pm 1 }$ . This gives us a very flexible way of constructing a wide variety of groups. We might, for example, use a group presentation to encode a question from elsewhere in mathematics. We could then ask about the properties of the group thus defined, and see what they had to tell us about our original problem.

# 3 Why Study Finitely Presented Groups?

Groups arise throughout the whole of mathematics as groups of automorphisms. These are maps from an object to itself that preserve all of the defining structure: two examples are the invertible linear maps [I.3 §4.2] from a vector space [I.3 §2.3] to itself, and the homeomorphisms from a topological space [III.90] to itself. Groups encapsulate the essence of symmetry and for this reason demand our attention. We are driven to understand their general nature, identify groups that deserve particular attention, and develop techniques for constructing new groups (from old ones, or from new ideas). And, reversing the process of abstraction, when given a group, we want to find concrete instances of it. For example, we might like to realize it as the group of automorphisms of some interesting object, with the aim of illuminating the nature of both the object and the group. (See the article on representation theory [IV.9] for more on this theme.)

# 3.1 Why Present Groups in Terms of Generators and Relations?

The short answer is that this is the form in which groups often “appear in nature.” This is particularly true in topology. Before looking at a general result that illustrates this point, let us examine a simple example. Consider the group D of all isometries of R that are generated by the reflections at the points 0, 1, and 2: that is, the group generated by the three functions $\alpha _ { 0 }$ , $\alpha _ { 1 }$ , and $\alpha _ { 2 } ,$ , which take x $\mathbf { t o } - x , 2 - x$ , and $_ { 4 - x , }$ , respectively. You may recognize this group to be the infinite dihedral group, and you may notice that the generator $\alpha _ { 2 }$ is superfluous, since it can be generated from $\alpha _ { 0 }$ and $\alpha _ { 1 }$ . But let us close our eyes to these observations as we let a presentation emerge from the action.

To this end, we choose an open interval U with the property that the images of U under the maps in D cover the whole of the real line, say $U = ( - \frac { 1 } { 2 } , \bar { \frac { 3 } { 2 } } )$ . Now let us record two pieces of data: the only elements of D (apart from the identity) that fail to move U completely off itself are $\alpha _ { 0 }$ and $\alpha _ { 1 } .$ and, among all products of length at most 3 in those two letters, the only nontrivial ones that act as the identity on R are $\alpha _ { 0 } ^ { 2 }$ and $\alpha _ { 1 } ^ { 2 } .$ . You may like to prove that $\smash { \langle \alpha _ { 0 } , \alpha _ { 1 } \ | \ \alpha _ { 0 } ^ { 2 } , \alpha _ { 1 } ^ { 2 } \rangle }$ is a presentation of D.

This is in fact a special case of a general result, which we now state. (The proof of it is somewhat involved.) Let X be a topological space that is both path connected [IV.6 §1] and simply connected [III.93], and let Γ be a group of homeomorphisms from X to itself. Then any choice of path-connected open subset $U \subset X$ such that the images of U cover all of X gives rise to a presentation ${ \cal { T } } = \langle S \mid R \rangle$ , where $S = \{ \gamma \in T \mid \gamma ( U ) \cap U \neq \emptyset \}$ and R consists of all words w  F(S) of length at most 3 such that $w = 1$ in Γ . Thus, the identification of a suitable subset U provides one with a presentation of Γ , and the task of a group theorist is to determine the nature of the group from this information.

To see how difficult this task is, you might like to consider the groups

$$
G _ {n} = \langle a _ {1}, \dots , a _ {n} \mid a _ {i} ^ {- 1} a _ {i + 1} a _ {i} a _ {i + 1} ^ {- 2}, i = 1, \dots , n \rangle ,
$$

where we interpret i 1 as 1 when i  n. One of $G _ { 3 }$ and $G _ { 4 }$ is trivial and the other is infinite. Can you decide which is which?

To illustrate a more subtle point, let us consider a finitely presented group that we perhaps feel we understand: the group $T _ { \Delta }$ that we were discussing earlier. If we want to describe this group to a blind friend unfamiliar with the triangular tiling of the plane, what can we say to make her understand the group, or at least convince her that we understand the group?

Our friend might reasonably ask us to list the elements of our group, so we begin to describe them as products (words) in the given generators. But as we begin to do so we hit a problem: we do not want to list any element more than once and in order to avoid redundancy we have to know which pairs of words $w _ { 1 }$ , w2 represent the same element of $I _ { \Delta } ;$ equivalently, we must be able to recognize which words $w _ { 1 } ^ { - 1 } w _ { 2 }$ are relations in the group. Determining which words are relations is called the word problem for the group. Even in ${ \cal { T } } _ { \Delta }$ this takes some work, and in the groups $G _ { n }$ we quickly find ourselves at a loss.

Note that as well as allowing one to list the elements of the group effectively, a solution to the word problem also allows one to determine the multiplication table, since deciding whether $w _ { 1 } w _ { 2 } = w _ { 3 }$ is the same as deciding whether $w _ { 1 } w _ { 2 } w _ { 3 } ^ { - 1 } = 1$ .

# 3.2 Why Finitely Presented Groups?

The packaging of infinite objects into finite amounts of data arises throughout mathematics in the various guises of compactness [III.9]. Finite presentation is basically a compactness condition: a group can be finitely presented if and only if it is the fundamental group of a reasonable compact space, as we shall see later.

Another good reason for studying finitely presented groups is that the Higman embedding theorem (to be discussed later) allows us to encode questions about arbitrary turing machines [IV.20 §1.1] as questions about such groups and their subgroups.

# 4 The Fundamental Decision Problems

In exploring the geometry and topology of low-dimensional manifolds at the beginning of the twentieth century, Max Dehn saw that many of the problems that he was wrestling with could be “reduced” to questions about finitely presented groups. For example, he gave a simple formula for associating with a knot diagram [III.44] a finite presentation of a group. There was one relation for each crossing in the diagram and he argued that the resulting group would be isomorphic to Z if and only if the knot was the unknot: that is, if and only if it could be continuously deformed into a circle. It is extremely hard to tell by staring at a knot diagram whether it is actually the unknot, so this seems like a useful reduction until one realizes that it can be just as hard to tell whether a finitely presented group is isomorphic to Z. For example, here is the presentation of Z that Dehn’s recipe associates with one of smallest possible pictures of the unknot, namely a diagram with just four crossings:

$$
\langle a _ {1}, a _ {2}, a _ {3}, a _ {4}, a _ {5} \mid
$$

$$
a _ {1} ^ {- 1} a _ {3} a _ {4} ^ {- 1}, a _ {2} a _ {3} ^ {- 1} a _ {1}, a _ {3} a _ {4} ^ {- 1} a _ {2} ^ {- 1}, a _ {4} a _ {5} ^ {- 1} a _ {4} a _ {3} ^ {- 1} \rangle .
$$

Thus Dehn’s investigations led him to understand how difficult it is to extract information from a group presentation. In particular, he was the first to identify the fundamental role of the word problem, which we alluded to earlier, and he was one of the first to begin to understand that there are fundamental problems associated with the challenge of developing algorithms that extract knowledge from well-defined objects such as group presentations. In his famous article of 1912 Dehn writes:

The general discontinuous group is given by n generators and m relations between them. … Here there are above all three fundamental problems whose solution is very difficult and which will not be possible without a penetrating study of the subject.

1. The identity [word] problem: An element of the group is given as a product of generators. One is required to give a method whereby it may be decided in a finite number of steps whether this element is the identity or not.   
2. The transformation [conjugacy] problem: Any two elements S and T of the group are given. A method is sought for deciding the question whether S and T can be transformed into each other, i.e., whether there is an element U of the group satisfying the relation

$$
S = U T U ^ {- 1}.
$$

3. The isomorphism problem: Given two groups, one is to decide whether they are isomorphic or not (and further, whether a given correspondence between the generators of one group and elements of the other is an isomorphism or not).

We shall take these problems as the starting point for three lines of enquiry. First, we shall work toward an outline of the proof that all of these problems are, in a strict sense, unsolvable for general finitely presented groups.

The second use that we shall make of Dehn’s problems is to hold them up as fundamental measures of complexity for each of the classes of groups that we subsequently encounter. If we can prove, for example, that the isomorphism problem is solvable in one class of groups but not in another, then we will have given genuine substance to previously vague assertions to the effect that the second class is “harder.”

Finally, I want to make the point that geometry lies at the heart of the fundamental issues in combinatorial group theory: it may not be immediately obvious, but its implicit presence is nonetheless a fundamental trait of group theory and not something imposed for reasons of taste. To illustrate this point I shall explain how the study of the large-scale geometry of leastarea disks in riemannian manifolds [I.3 §6.10] is intimately connected with the study of the complexity of word problems in arbitrary finitely presented groups.

# 5 New Groups from Old

Suppose that you have two groups, $G _ { 1 }$ and $G _ { 2 }$ , and want to combine them to form a new group. The first method that is taught in a typical course on group theory is to take the Cartesian product $G _ { 1 } \times G _ { 2 } \colon$ a typical element has the form $( g , h )$ with $g ~ \in ~ G _ { 1 }$ and $h \in G _ { 2 }$ , and the product of $( g , h )$ with $( g ^ { \prime } , h ^ { \prime } )$ is defined to be $( g g ^ { \prime } , h h ^ { \prime } )$ . The set of elements of the form $( g , e )$ (where e is the identity of $G _ { 2 } )$ is a copy of $G _ { 1 }$ inside $G _ { 1 } \times G _ { 2 } ,$ and similarly the set of elements of the form $( e , h )$ is a copy of $G _ { 2 }$ .

These copies have nontrivial relations between their elements: for example, $( e , h ) ( g , e ) \ = \ ( g , e ) ( e , h )$ . We would now like to take two groups $T _ { 1 }$ and $T _ { 2 }$ and combine them in a different way to form a group called the free product $I _ { 1 } * I _ { 2 }$ , which contains copies of ${ \cal { T } } _ { 1 }$ and $T _ { 2 }$ and as few additional relations as possible. That ${ \mathrm { i } } s ,$ we would like there to be embeddings $i _ { j } : { \cal T } _ { j } \hookrightarrow { \cal T } _ { 1 } \ast { \cal T } _ { 2 }$ so that $i _ { 1 } ( I _ { 1 } )$ and $i _ { 2 } ( I _ { 2 } )$ generate $I _ { 1 } * I _ { 2 }$ but they are not intertwined in any way. This requirement is neatly encapsulated by the following universal property: given any group G and any two homomorphisms $\phi _ { 1 } : I _ { 1 } \to G$ and $\phi _ { 2 } : I _ { 2 } \to G$ , there should be a unique homomorphism $\phi : I _ { 1 }$ ∗ $\Gamma _ { 2 }  G$ such that $\Phi \circ i _ { j } = \phi _ { j } \mathrm { f o r } \ j = 1 , 2$ . (Less formally, Φ behaves like $\phi _ { 1 }$ on the copy of $\boldsymbol { { \cal T } _ { 1 } }$ and behaves like φ2 on the copy of $I _ { 2 } . )$

It is easy to check that this property characterizes ${ { T } _ { 1 } } ~ * ~ { { T } _ { 2 } }$ up to isomorphism, but it leaves open the question of whether $I _ { 1 } * I _ { 2 }$ actually exists. (These are the standard pros and cons of defining an object by means of a universal property.) In the present setting, existence is easily established using presentations: let $\langle A _ { 1 } \ | \ R _ { 1 } \rangle$ be a presentation of $\boldsymbol { { \cal T } _ { 1 } }$ and let $\langle A _ { 2 } \ | \ R _ { 2 } \rangle$ be a presentation of $\displaystyle { \cal { I } } _ { 2 } ,$ , with $A _ { 1 }$ and $A _ { 2 }$ disjoint, and then define $I _ { 1 } * I _ { 2 }$ to be $\langle A _ { 1 } \sqcup A _ { 2 } \ | \ R _ { 1 } \sqcup R _ { 2 } \rangle$ (where % denotes a union of disjoint sets).

More intuitively, one can define ${ \cal { T } } _ { 1 }$ \* $T _ { 2 }$ to be the set of alternating sequences $a _ { 1 } b _ { 1 } \cdots a _ { n } b _ { n }$ with each $a _ { i }$ belonging to $\boldsymbol { { \cal T } _ { 1 } }$ and each $b _ { j }$ belonging to $T _ { 2 }$ , with the extra condition that none of the $\alpha _ { i }$ and $b _ { j }$ equals the identity, except possibly $_ { a _ { 1 } }$ or $b _ { n }$ . The group operations in ${ \cal { T } } _ { 1 }$ and $T _ { 2 }$ extend to this set in an obvious way: for example, $( a _ { 1 } b _ { 1 } a _ { 2 } ) ( a _ { 1 } ^ { \prime } b _ { 1 } ^ { \prime } ) = a _ { 1 } b _ { 1 } a _ { 2 } ^ { \prime } b _ { 1 } ^ { \prime }$ , where $a _ { 2 } ^ { \prime } = a _ { 2 } a _ { 1 } ^ { \prime }$ , except that if $a _ { 2 } a _ { 1 } ^ { \prime } = 1$ then the product cancels down to $a _ { 1 } b _ { 2 } ^ { \prime }$ , where $b _ { 2 } ^ { \prime } = b _ { 1 } b _ { 1 } ^ { \prime }$ .

Free products occur naturally in topology: if one has topological spaces $X _ { 1 } , X _ { 2 }$ with marked points $p _ { 1 } \in X _ { 1 } ,$ $p _ { 2 } \in X _ { 2 }$ , then the fundamental group [IV.6 §2] of the space $X _ { 1 } \lor X _ { 2 }$ obtained from $X _ { 1 } \sqcup X _ { 2 }$ by making the identification $p _ { 1 } = p _ { 2 }$ is the free product of $\pi _ { 1 } ( X _ { 1 } , p _ { 1 } )$ ) and $\pi _ { 1 } ( X _ { 2 } , p _ { 2 } )$ . The Seifert–van Kampen theorem tells one how to present the fundamental group of a space obtained by gluing $X _ { 1 }$ and $X _ { 2 }$ along larger subspaces. If the inclusion of the subspaces gives rise to an injection of fundamental groups, then one can express the fundamental group of the resulting space as an amalgamated free product, which we now define.

Let ${ \cal { T } } _ { 1 }$ and $T _ { 2 }$ be two groups. If some other group contains copies of $T _ { 1 }$ and $T _ { 2 }$ , then the intersection of those copies must contain the identity element. The free product ${ { T } _ { 1 } } * { { T } _ { 2 } }$ was the freest group we could build that was subject to this minimal constraint. Now we shall insist that the copies of ${ \cal { T } } _ { 1 }$ and $T _ { 2 }$ intersect nontrivially, specify which of their subgroups must lie in the intersection, and build the freest group that satisfies this constraint.

Suppose, then, that $A _ { 1 }$ is a subgroup of ${ \cal { T } } _ { 1 }$ and that φ is an isomorphism from $A _ { 1 }$ to a subgroup $A _ { 2 }$ of $T _ { 2 }$ . As in the example of the free product, one can define the “freest product that identifies $A _ { 1 }$ and $A _ { 2 } ^ { \prime \prime }$ by means of a universal property. Again, one can establish the existence of such a group using presentations: if ${ \cal { I } } _ { 1 } =$ $\langle S _ { 1 } \mid R _ { 1 } \rangle$ and ${ \cal T } _ { 2 } = \langle S _ { 2 } \mid { \cal R } _ { 2 } \rangle$ , the group we seek takes the form

$$
\langle S _ {1} \sqcup S _ {2} \mid R _ {1} \sqcup R _ {2} \sqcup T \rangle .
$$

Here, $T = \{ u _ { a } \nu _ { a } ^ { - 1 } \mid a \in A _ { 1 } \}$ , where $u _ { a }$ is some word that represents a in (the presentation of) $\boldsymbol { { \cal T } _ { 1 } }$ and $\nu _ { a }$ is a word that represents $\phi ( a )$ in $T _ { 2 }$ .

This group is called the amalgamated free product of $\boldsymbol { { \cal T } _ { 1 } }$ and $T _ { 2 }$ along $A _ { 1 }$ and $A _ { 2 }$ . It is often described by the casual and ambiguous notation $\boldsymbol { { \cal T } } _ { 1 } * _ { A _ { 1 } = A _ { 2 } } \boldsymbol { { \cal T } } _ { 2 }$ , or even $\Gamma _ { 1 } * _ { A } { \Gamma } _ { 2 }$ , where $A \cong A _ { j }$ is an abstract group.

Unlike with free products, it is no longer obvious that the maps ${ { T } _ { i } }  { { T } _ { 1 } } * _ { A } { { T } _ { 2 } }$ implicit in this construction are injective, but they do turn out to be, as was shown by Schreier in 1927.

A related construction of Higman, Neumann, and Neumann in 1949 answers the following question: given a group $\boldsymbol { { \cal T } }$ and an isomorphism $\psi : B _ { 1 } \ \to \ B _ { 2 }$ between subgroups of $T ,$ can one always embed Γ in a bigger group so that ψ becomes the restriction to $B _ { 1 }$ of a conjugation?

By now, having seen the idea in the context of both free products and amalgamated free products, the reader may guess how one goes about answering this question: one writes down the presentation of a universal candidate for the desired enveloping group, denoted $T * _ { \psi } { } _ { ; }$ , and then one sets about proving that the natural map from $\boldsymbol { { \cal T } }$ to $T * _ { \psi }$ (which takes each word to itself) is injective. Thus, given ${ \cal { T } } = \langle A \mid { \cal { R } } \rangle$ , we introduce a symbol t ∉ A (usually called the stable letter), we choose for each $b \in B _ { 1 }$ words $\hat { b } , \tilde { b } \in F ( A )$ with ${ \hat { b } } = b$ and ${ \tilde { b } } = \psi ( b )$ in Γ , and we define

$$
\Gamma * _ {\psi} = \langle A, t \mid R, t \hat {b} t ^ {- 1} \tilde {b} ^ {- 1} (b \in B _ {1}) \rangle .
$$

This is the freest group we can build from Γ by adjoining a new element t and requiring it to satisfy all the equations we want it to, namely $t \hat { b } t ^ { - 1 } = \tilde { b }$ for every $b \in B _ { 1 }$ (which we can think of as saying that $t b t ^ { - 1 } =$ $\psi ( b ) )$ . This group is called an HNN extension of Γ (after Higman, Neumann, and Neumann).

Now we must show that the natural map from Γ to $T * _ { \psi }$ is injective. That is, if you take an element γ of Γ and regard it as an element of $T * _ { \psi } ,$ , you should not be able to use t and the relations in $T * _ { \psi }$ to cancel γ down to the identity. This is proved with the help of the following more general result known as Britton’s lemma. Suppose that w is a word in the free group $F ( A , t )$ . Then the only circumstances under which it can give rise to the identity in the group $T * _ { \psi }$ are if either it does not involve t and represents the identity in Γ or it involves t but can be simplified in an obvious way by containing a “pinch.” A pinch is a subword of the form $t b t ^ { - 1 }$ , where b is a word in $F ( A )$ ) that represents an element of $B _ { 1 }$ (in which case we can replace it by $\psi ( b ) )$ , or one of the form $t ^ { - 1 } b ^ { \prime } t$ , where $b ^ { \prime }$ represents an element of $B _ { 2 }$ (in which case we can replace it by $\psi ^ { - 1 } ( b ^ { \prime } ) )$ . Thus, if you are given a word that involves t and contains no pinches, then you know that it cannot be canceled down to the identity.

A similar noncancellation result holds for the amalgamated free product $T _ { 1 } * _ { A _ { 1 } = A _ { 2 } } * _ { 2 } . \mathrm { I f } \ g _ { 1 } , \dots , g _ { n }$ belong to $T _ { 1 }$ but not to $A _ { 1 }$ and $h _ { 1 } , \ldots , h _ { n }$ belong to $T _ { 2 }$ but not to $A _ { 2 } ,$ then the word $g _ { 1 } h _ { 1 } g _ { 2 } h _ { 2 } \cdot \cdot \cdot \cdot g _ { n } h _ { n }$ cannot equal the identity in $\boldsymbol { { \cal T } } _ { 1 } * _ { A _ { 1 } = A _ { 2 } } \boldsymbol { { \cal T } } _ { 2 }$ .

These noncancellation results do far more than show that the natural homomorphisms we have been considering are injective: they also demonstrate further aspects of freeness in amalgamated free products and HNN extensions. For example, suppose that in the amalgamated free product $\boldsymbol { { \cal T } } _ { 1 } * _ { A _ { 1 } = A _ { 2 } } \boldsymbol { { \cal T } } _ { 2 }$ we can find an element $_ g$ of $T _ { 1 }$ that generates an infinite group that intersects $A _ { 1 }$ in the identity and an element h of $T _ { 2 }$ that does the same for $A _ { 2 }$ . Then the subgroup of $\begin{array} { r } { { \cal T } _ { 1 } * _ { A _ { 1 } = A _ { 2 } } { \cal T } _ { 2 } } \end{array}$ generated by $_ g$ and h is the free group on those two generators. With a little more effort, one can deduce that any finite subgroup of $\boldsymbol { { \cal T } } _ { 1 } * _ { A _ { 1 } = A _ { 2 } } \boldsymbol { { \cal T } } _ { 2 }$ has to be conjugate to a subgroup of the obvious copy of either $T _ { 1 }$ or $\displaystyle { \cal { I } } _ { 2 } .$ . Similarly, the finite subgroups of $T * _ { \psi }$ are conjugates of subgroups of Γ . We shall exploit these facts in the constructions that follow.

There are many ways of combining groups that I have not mentioned here. I have chosen to focus on amalgamated free products and HNN extensions partly because they lead to transparent solutions of the basic problems discussed below but more because of their primitive appeal and the way in which they arise naturally in the calculation of fundamental groups. They also mark the beginning of arboreal group theory, which we will discuss later. If space allowed, I would ${ \bf g 0 }$ on to describe semidirect and wreath products, which are also indispensable tools of the group theorist.

Before turning to some applications of HNN extensions and amalgamated free products, I want to return to the Burnside problem, which asks if there exist finitely generated infinite groups all of whose elements have a given finite order. This question generated important developments throughout the twentieth century, particularly in Russia. It is appropriate to mention it here because it provides another illustration of the fact that it can be useful to study a universal object in order to solve a general question.

# 5.1 The Burnside Problem

Given an exponent m, one clarifies the problem at hand by considering the free Burnside group $B _ { n , m }$ given by the presentation $\langle a _ { 1 } , \dots , a _ { n } \mid R _ { m } \rangle$ , where $R _ { m }$ consists of all mth powers in the free group $F ( a _ { 1 } , \ldots , a _ { n } )$ . It is clear that $B _ { n , m }$ maps onto any group with at most n generators in which every element has order dividing m. Therefore, there exists a finitely generated infinite group with all elements of the same finite order if and only if, for suitable values of n and $m$ , the group $B _ { n , m }$ is infinite. Thus, a question that takes the form, Does there exist a group such that $\cdots 3 ,$ becomes a question about just one group.

Novikov and Adian showed in 1968 that $B _ { n , m }$ is infinite when $n \geqslant 2$ and m $\geqslant ~ 6 6 7$ is odd. Determining the exact range of values for which $B _ { n , m }$ is infinite is an active area of research. Of far greater interest is the open question of whether there exist finitely presented infinite groups that are quotients of $B _ { n , m }$ . Zelmanov was awarded the Fields Medal for proving that each $B _ { n , m }$ has only finitely many finite quotients.

# 5.2 Every Countable Group Can Be Embedded in a Finitely Generated Group

Given a countable group G we can list its elements, $g _ { 0 } , g _ { 1 } , g _ { 2 } , \ldots$ , taking $g _ { 0 }$ to be the identity. We can then take a free product of G with an infinite cyclic group $\langle s \rangle \cong \mathbb { Z } . \operatorname { L e t } \Sigma _ { 1 }$ be the set of all elements of $G * \mathbb { Z }$ of the form $s _ { n } = g _ { n } s ^ { n }$ with $n \geqslant 1$ . Then the subgroup $\langle \Sigma _ { 1 } \rangle$ generated by $\Sigma _ { 1 }$ is isomorphic to the free group $F ( \Sigma _ { 1 } )$ . Similarly, if we let $\Sigma _ { 2 } = \{ s _ { 2 } , s _ { 3 } , \dots \}$ (so it is $\Sigma _ { 1 }$ with the element $s _ { 1 } = g _ { 1 } s$ removed), then $\langle \Sigma _ { 2 } \rangle$ is isomorphic to $F ( \Sigma _ { 2 } )$ . It follows that the map $\psi ( s _ { n } ) = s _ { n + 1 }$ gives rise to an isomorphism from $\langle \Sigma _ { 1 } \rangle$ to $\left. \Sigma _ { 2 } \right.$ . Now take the HNN extension $( G * \mathbb { Z } ) * _ { \psi }$ , whose stable letter we denote by t. This group contains a copy of $G ,$ as we noted before. Moreover, since we have ensured that $t s _ { n } t ^ { - 1 } = s _ { n + 1 }$ for every $n \geqslant 1$ , it can be generated by just the three elements $s _ { 1 } , s ,$ , and t. Thus, we have embedded an arbitrary countable group into a group with three generators. (We leave the reader to think about how one can vary this construction to produce a group with two generators.)

# 5.3 There Are Uncountably Many Nonisomorphic Finitely Generated Groups

This was proved by B. H. Neumann in 1932. Since there are infinitely many primes, there are uncountably many nonisomorphic groups of the form $\oplus _ { p \in P } \mathbb { Z } _ { p }$ , where $P$ is an infinite set of primes. We have seen that each of these groups can be embedded in a finitely generated group, and our earlier comments on finite subgroups of HNN extensions show that no two of the resulting finitely generated groups are isomorphic.

# 5.4 An Answer to Hopf’s Question

A group G is called Hopfian if every surjective homomorphism from G to G is an isomorphism. Most familiar groups have this property: for example, finite groups obviously do, as do $\mathbb { Z } ^ { n }$ (as you can prove using linear algebra) and free groups. So too do groups of matrices such as ${ \mathrm { S L } } _ { n } ( \mathbb { Z } )$ , as we shall discuss in a moment. A simple example of a non-Hopfian group is the group consisting of all infinite sequences of integers (under pointwise addition), since the function that takes $( a _ { 1 } , a _ { 2 } , a _ { 3 } , \dots )$ to $( a _ { 2 } , a _ { 3 } , a _ { 4 } , \dots )$ is a surjective homomorphism that contains $( 1 , 0 , 0 , \ldots )$ in its kernel. But is there a finitely presented example? The answer is yes, and Higman was the first to construct one. The following examples are due to Baumslag and Solitar.

Let $p \geqslant$ 2 be an integer and identify Z with the free group a generated by a single generator a. Then the subgroups $p \mathbb { Z }$ and $( p + 1 ) \mathbb { Z } { \mathbf o } { \mathbf f } \ \mathbb { Z }$ are identified with the powers of $\boldsymbol { a } ^ { p }$ and $a ^ { p + 1 }$ , respectively. Let $\psi$ be the isomorphism between these subgroups that takes $a ^ { p }$ to $a ^ { p + 1 }$ and consider the corresponding HNN extension

B. This has presentation $B = \langle a , t \mid t a ^ { - p } t ^ { - 1 } a ^ { p + 1 } \rangle$ . The homomorphism $\psi : B  B$ defined by $t \mapsto t , a \mapsto a ^ { p }$ is clearly a surjection but its kernel contains, for example, the element $c = a t a ^ { - 1 } t ^ { - 1 } a ^ { - 2 } t a t ^ { - 1 } a ,$ , which does not contain a pinch and is therefore not equal to the identity, by Britton’s lemma. (If you want to convince yourself how useful this lemma is, set $p = 3$ and try to prove directly that c is not equal to the identity in the group B just defined.)

# 5.5 A Group That Has No Faithful Linear Representation

One can show that a finitely generated group G of matrices over any field is residually finite, which means that for each nontrivial element $g \in G$ there exists a finite group Q and a homomorphism $\pi : G \to Q$ with $\pi ( g ) \neq 1$ . For example, if you are given an element $g \in { \mathrm { S L } } _ { n } ( \mathbb { Z } )$ , then you can pick an integer m bigger than the absolute values of all the entries in $_ g$ (which is an n  n matrix) and consider the homomorphism from ${ \mathrm { S L } } _ { n } ( \mathbb { Z } )$ to $\mathrm { S L } _ { n } ( \mathbb { Z } / m \mathbb { Z } )$ that reduces the matrix entries mod m. The image of g in the finite group $\mathrm { S L } _ { n } \left( \mathbb { Z } / m \mathbb { Z } \right)$ is clearly nontrivial.

Non-Hopfian groups are not residually finite, and hence are not isomorphic to a group of matrices over any field. One can see that the non-Hopfian group B defined above is not residually finite by considering what happens to the nontrivial element c. We saw that there was a surjective homomorphism $\psi : B  B$ with $\psi ( c ) = 1$ . Let $c _ { n }$ be an element such that $\psi ^ { n } ( c _ { n } ) = c$ (which exists since ψ is a surjection). If there were a homomorphism π from B to a finite group Q with $\pi ( c ) \neq 1$ , then we would have infinitely many distinct homomorphisms from B to Q, namely the compositions $\pi \circ \psi ^ { n }$ ; these are distinct because π ◦ $\psi ^ { m } ( c _ { n } ) = 1$ if m > n and $\pi \circ \psi ^ { n } ( c _ { n } ) = \pi ( c ) \neq 1$ . This is a contradiction, since a homomorphism from a finitely generated group to a finite group is determined by what it does to the generators, so there can only be finitely many such homomorphisms.

# 5.6 Infinite Simple Groups

Britton’s lemma actually tells us more than that $c \neq 1$ : the subgroup Λ of B generated by t and c is in fact a free group on those generators. Thus we may form the amalgamated free product Γ of two copies of $B ,$ denoted $B _ { 1 }$ and $B _ { 2 }$ , by gluing together the two copies of Λ with the isomorphism $c _ { 1 } \mapsto t _ { 2 } , t _ { 1 } \mapsto c _ { 2 }$ . We have seen that in any finite quotient of $I = B _ { 1 } * _ { \varLambda } B _ { 2 }$ , the elements $c _ { 1 } ~ ( = t _ { 2 } )$ and $c _ { 2 } ~ ( = t _ { 1 } )$ must have trivial image, and it is easy to deduce from this that in fact the quotient must be trivial. Thus Γ is an infinite group with no finite quotients. It follows that the quotient of Γ by any maximal proper normal subgroup is also infinite (and it is simple by maximality).

The simple group that we have constructed is infinite and finitely generated but it is not finitely presentable. Finitely presented infinite simple groups do exist, but they are much harder to construct.

# 6 Higman’s Theorem and Undecidability

We have seen that there are uncountably many (nonisomorphic) finitely generated groups. But as there are only countably many finitely presented groups, only countably many finitely generated groups can be subgroups of finitely presented groups. Which ones are they?

A complete answer to this question is provided by a beautiful and deep theorem proved by Graham Higman in 1961, which says, roughly, that the groups that arise are all those that are algorithmically describable. (If you have no idea what this means, even roughly, then you might like to read the insolubility of the halting problem [V.20] before continuing with this section.)

A set S of words over a finite alphabet A is called recursively enumerable if there is some algorithm (or more formally, Turing machine) that can produce a complete list of the elements of S. A case of particular interest is when A is just a singleton, in which case a word is determined by its length and we can think of S as a set of nonnegative integers. The elements of S need not be listed in a sensible order, so having an algorithm that produces an exhaustive list of S does not mean that one can use the algorithm to determine that some given word w does not belong to S: if you imagine standing by your computer as it enumerates S, there will not in general come a time when you can say to yourself, “If it was going to appear, then it would have done so by now,” and therefore be certain that it is not in S. If you want an algorithm with this further property, then you need the stronger notion of a recursive set, which is a set S such that S and its complement are both recursively enumerable. Then you can list all the elements that belong to S and you can also list all the elements that do not belong to S.

A finitely generated group is said to be recursively presentable if it has a presentation with a finite number of generators and a recursively enumerable set of defining relations. In other words, such a group is not necessarily finitely presented, but at least the presentation of the group is “nice” in the sense that it can be generated by some algorithm.

Higman’s embedding theorem states that a finitely generated group G is recursively presentable if and only if it is isomorphic to a subgroup of a finitely presented group.

To get a feeling for how nonobvious this is, you might consider the following presentation of the group of all rationals under addition, in which the generator $a _ { n }$ corresponds to the fraction 1/n!:

$$
Q = \left\langle a _ {1}, a _ {2}, \dots \mid a _ {n} ^ {n} = a _ {n - 1} \forall n \geqslant 2 \right\rangle .
$$

Higman’s theorem tells us that Q can be embedded in a finitely presented group, but no truly explicit embedding is known.

The power of Higman’s theorem is illustrated by the ease with which it implies the celebrated undecidability results that were rightly regarded as watersheds of twentieth-century mathematics. In order to make this case convincingly, I shall give a complete proof (except that I shall assume some of the facts mentioned earlier) that there exist finitely presented groups with unsolvable word problems, and also that there are sequences of finitely presented groups among which one cannot decide isomorphism. We shall also see how these group-theoretic results can be used to translate undecidability phenomena into topology.

The basic seed of undecidability comes from the fact that there are recursively enumerable subsets $S \subset \mathbb { N }$ that are not recursive. Using this fact one can readily construct finitely generated groups with an unsolvable word problem: given such a set of integers S we consider

$$
J = \langle a, b, t \mid t (b ^ {n} a b ^ {- n}) t ^ {- 1} = b ^ {n} a b ^ {- n} \forall n \in S \rangle .
$$

This is the HNN extension of the free group $F ( a , b )$ associated with the identity map $L \to L ,$ , where L is the subgroup generated by $\{ b ^ { n } a b ^ { - n } : n \in S \}$ . Britton’s lemma tells us that the word

$$
w _ {m} = t (b ^ {m} a b ^ {- m}) t ^ {- 1} (b ^ {m} a ^ {- 1} b ^ {- m})
$$

equals $1 \in J$ if and only if $m \in S ,$ and by definition there is no algorithm to decide if m S, so we cannot decide which of the $\boldsymbol { w } _ { m }$ are relations. Thus J has an unsolvable word problem.

That there exist finitely presented groups for which the word problem is unsolvable is a much deeper fact, but with Higman’s embedding theorem at hand the proof becomes almost trivial: Higman tells us that J can be embedded in a finitely presented group Γ , and it is a relatively straightforward exercise to show that if one cannot decide which words in the generators of J represent the identity, then one cannot decide for arbitrary words in the generators of Γ either.

Once one has a finitely presented group with an unsolvable word problem, it is easy to translate undecidability into all manner of other problems. For example, suppose that ${ \cal { T } } ~ = ~ \langle { A } ~ | ~ R \rangle$ is a finitely presented group with an unsolvable word problem, where A = $\{ a _ { 1 } , \ldots , a _ { n } \}$ and no ai equals the identity in Γ . For each word w made out of the letters in A and their inverses, define a group $\varGamma _ { w }$ to have presentation

$$
\langle A, s, t \mid R, t ^ {- 1} (s ^ {i} a _ {i} s ^ {- i}) t (s ^ {i} w s ^ {- i}), i = 1, \dots , n \rangle .
$$

It is not hard to show that if $w = 1$ in Γ then $\varGamma _ { w }$ is the free group generated by s and t. If $w \ne 1$ , then $\varGamma _ { w }$ is an HNN extension. In particular, it contains a copy of Γ , and hence has an unsolvable word problem, which means that it cannot be a free group. Thus, since there is no algorithm to decide whether $w = 1$ in Γ , one cannot decide which of the groups $\varGamma _ { w }$ are isomorphic to which others.

A variant of this argument shows that there is no algorithm to determine whether or not a given finitely presented group is trivial.

We shall see in a moment that every finitely presented group G is the fundamental group of some compact four-dimensional manifold. By following a standard proof of this theorem with considerable care, Markov proved in 1958 that in dimensions 4 and above there is no algorithm to decide which compact manifolds (presented as simplicial complexes, for example) are homeomorphic. His basic idea was to show that if there were an algorithm to determine which triangulated 4-manifolds are homeomorphic, then one could use it to determine which finitely presented groups are trivial, which we know is impossible. In order to implement this idea one has to be careful to arrange that the 4-manifolds associated with different presentations of the trivial group are homeomorphic: this is the delicate part of the argument.

Strikingly, there does exist an algorithm to decide which compact three-dimensional manifolds are isomorphic. This is an extremely deep theorem that relies in particular on Perelman’s solution to thurston’s geometrization conjecture [IV.7 §2.4].

# 7 Topological Group Theory

Let us change perspective now and look at the symbols $\boldsymbol { P } \equiv \langle a _ { 1 } , \dots , a _ { 2 } \mid r _ { 1 } , \dots , r _ { m } \rangle$ through the eyes of a topologist. Instead of interpreting P as a recipe for constructing a group, we regard it as a recipe for constructing a topological space [III.90], or more specifically a two-dimensional complex. Such spaces consist of points, called vertices, some of which are linked by directed paths, called edges, or 1-cells. If a collection of such 1-cells forms a cycle, then it can be filled in with a face, or 2-cell: topologically speaking, each face is a disk with a directed cycle as its boundary.

To see what this complex is, let us first consider the standard presentation $P \ \equiv \ \langle a , b \mid a b a ^ { - 1 } b ^ { - 1 } \rangle$ of $ { \mathbb { Z } } ^ { 2 }$ . (This is generated by a and b and the relation tells us that $a b \ = \ b a . )$ We begin with a graph $K ^ { 1 }$ that has a single vertex and two edges (which are loops) that are directed and labeled a and b. Next, we take a square $[ 0 , 1 ] \times [ 0 , 1 \dot { }$ ], the sides of which are directed and labeled $a , b , a ^ { - 1 } , b ^ { - 1 }$ as we proceed around the boundary. Imagine gluing the boundary of the square to the graph so as to respect the labeling of edges: with a bit of thought, you should be able to see that the result is a torus, that is, a surface in the shape of a bagel. An observation that turns out to be important is that the fundamental group of the torus is $\mathbb { Z } ^ { 2 }$ , the group we started with.

The idea of “gluing” is made precise by the use of attaching maps: we take a continuous map φ from the boundary of the square S to the graph $K ^ { 1 }$ that sends the corners of the square to the vertex of $K ^ { 1 }$ and sends each side (minus its vertices) homeomorphically onto an open edge. The torus is then the quotient of $K ^ { 1 } \sqcup S$ by the equivalence relation that identifies each x in the boundary of the square with its image $\phi ( x )$ .

With this more abstract language in hand, it is easy to see how the above construction generalizes to arbitrary presentations: given a presentation $P \equiv \langle a _ { 1 } , \ldots , a _ { n } \ |$ $r _ { 1 } , \ldots , r _ { m } \rangle$ , one takes a graph with a single vertex and n oriented loops, which are labeled $a _ { 1 } , \ldots , a _ { n }$ . Then for each $r _ { j }$ one attaches a polygonal disk by gluing its boundary circuit to the sequence of oriented edges that traces out the word $r _ { j }$ .

In general, the result will not be a surface as it was for $\langle a , b \mid a b a ^ { - 1 } b ^ { - 1 } \rangle$ . Rather, it will be a two-dimensional complex with singularities along the edges and at the vertex. You may find it instructive to do some more examples. From $\langle \alpha \mid a ^ { 2 } \rangle$ one gets the projective plane; from $\langle a , b , c , d \mid a b a ^ { - 1 } b ^ { - 1 } , c d c ^ { - 1 } d \rangle$ one gets a torus and a Klein bottle stuck together at a point. Picturing the 2-complex for $\langle a , b \mid a ^ { 2 } , b ^ { 3 } , ( a b ) ^ { 3 } \rangle$ is already rather difficult.

The construction of $K ( P )$ is the beginning of topological group theory. The Seifert–van Kampen theorem (mentioned earlier) implies that the fundamental group of $K ( P )$ is the group presented by P. But the group no longer sits inertly in the form of an inscrutable presentation—now it acts on the universal covering [III.93] of $K ( P )$ by homeomorphisms known as “deck transformations.” Thus, through the simple construction of $K ( P )$ (and the elegant theory of covering spaces in topology) we achieve our aim of realizing an abstract finitely presented group as the group of symmetries of an object with a potentially rich structure, on which we can bring global geometric and topological techniques to bear.

To obtain an improved topological model for our group, we can embed $K ( P )$ in $\mathbb { R } ^ { 5 }$ (just as one can embed a finite graph [III.34] in R3) and consider the compact four-dimensional manifold M obtained by taking all points that are a small fixed distance from the image. (I am assuming that the embedding is suitably “tame,” which one can arrange.) The mental picture to strive for here is a higher-dimensional analogue of the surface (sleeve) one gets by taking the points in $\mathbb { R } ^ { 3 }$ that are a small fixed distance from an embedded graph. The fundamental group of M is again the group presented by P, so now we have our arbitrary finitely presented group acting on a manifold (the universal cover of M). This allows us to use the tools of analysis and differential geometry.

The constructions of K(P ) and M establish the more difficult implication of the theorem, promised earlier, that a group can be finitely presented if and only if it is the fundamental group of a compact cell complex and of a compact 4-manifold. This result raises several natural questions. First, are there better, more informative, topological models for an arbitrary finitely presented group Γ ? And if not, then what can one say about the classes of groups defined by the natural constraints that arise when one tries to improve the model? For example, we would like to construct a lower-dimensional manifold with fundamental group Γ , enabling us to exploit our physical insight into threedimensional geometry. But it turns out that the fundamental groups of compact three-dimensional manifolds are very special; this observation lies near the heart of a great deal of mathematics at the end of the twentieth century. Other interesting fields open up when one asks which groups arise as the fundamental groups of compact spaces satisfying curvature [III.13] conditions, or constraints coming from complex geometry.

A particularly rich set of constraints comes from the following question. Can one arrange for an arbitrary finitely presented group to be the fundamental group of a compact space (a complex or manifold, perhaps) whose universal cover is contractible [IV.6 §2]? This is a natural question from the point of view of topology because a space with a contractible universal cover is, up to homotopy [IV.6 §2], completely determined by its fundamental group. If the fundamental group is $T ,$ then such a space is called a classifying space for Γ and its homotopy-invariant properties provide a rich array of invariants for the group Γ (getting away from the gross dependence that K(P ) has on P rather than Γ ).

If our earlier discussion of how hard it is to recognize Γ from P has left you very skeptical about whether this dependence can actually be removed, then your skepticism is well-founded: there are many obstructions to the construction of compact classifying spaces for an arbitrary finitely presented group; the study of them (under the generic name finiteness conditions) is a rich area at the interface of modern group theory, topology, and homological algebra.

One aspect of this area is the search for natural conditions that ensure the existence of compact classifying spaces (not necessarily manifolds). This is one of several places where manifestations of nonpositive curvature play a fundamental role in modern group theory. More combinatorial conditions also arise. For example, Lyndon proved that for any presentation $P \equiv$ $\langle A \mid r \rangle$ where the single defining relation $r \in F ( A )$ is not a nontrivial power, the universal cover of K(P ) is contractible.

A neighboring and highly active area of research concerns questions of uniqueness and rigidity for classifying spaces. (Here, as is common, the word rigidity is used to describe a situation in which requiring two objects to be equivalent in an apparently weak sense forces them to be equivalent in an apparently stronger sense.) For example, the (open) Borel conjecture asserts that if two compact manifolds have isomorphic fundamental groups and contractible universal covers, then those manifolds must be homeomorphic.

I have been talking mostly about realizing groups as fundamental groups, which led to certain free actions. That is, we could interpret the elements of the group as symmetries of a topological space and none of these symmetries had any fixed points. Before moving on to geometric group theory I should point out that there are many situations in which the most illuminating actions of a group are not free: one instead allows wellunderstood stabilizers. (The stabilizer of a point is the set of all symmetries in the group that leave that point fixed.) For example, the natural way in which to study ${ \cal { T } } _ { \Delta }$ is by its action on the triangulated plane, each vertex of which is left unmoved by twelve symmetries.

A deeper illustration of the merits of seeking insight into algebraic structure through nonfree actions on suitable topological spaces comes from the Bass–Serre theory of groups acting on trees, which subsumes the theory of amalgamated free products and HNN extensions, whose potency we saw earlier. (This theory and its extensions often go under the heading of arboreal group theory.)

A tree is a connected graph that has no circuits in it. It is helpful to regard it as a metric space [III.56] in which each edge has length 1. The group actions that one allows on trees are those that take edges to edges isometrically, never flipping an edge.

If a group Γ acts on a set X (in other words, if it can be regarded as a group of symmetries of X), then the orbit of a point $x \in X$ is the set of all its images gx with $g \in T . \mathrm { A }$ group Γ can be expressed as an amalgamated free product $A * _ { C } B \operatorname* { i f }$ and only if it acts on a tree in such a way that there are two orbits of vertices, one orbit of edges, and stabilizers A, B, C (where A and B are the stabilizers of adjacent vertices and intersect in C, which is the edge stabilizer). HNN extensions correspond to actions with one orbit of vertices and one orbit of edges. Thus, amalgamated free products and HNN extensions appear as graphs of groups, which are the basic objects of Bass–Serre theory. These objects allow one to recover groups acting on trees from the quotient data of the action, i.e., the quotient space (which is a graph) and the pattern of edge and vertex stabilizers.

An early benefit of Bass–Serre theory is a transparent and instructive proof that any finite subgroup of A C B is conjugate to a subgroup of either A or B: given any set V of vertices in a tree, there is a unique vertex or midpoint x minimizing max $\{ d ( x , \nu ) \mid \nu \in V \}$ ; one applies this observation with V an orbit of the finite subgroup; x provides a fixed point for the action of the subgroup; and any point stabilizer is conjugate to a subgroup of either A or B.

Arboreal group theory goes much deeper than this first application suggests. It is the basis for a decomposition theory of finitely presented groups from which it emerges, for example, that there is an essentially canonical maximal splitting of an arbitrary finitely presented group as a graph of groups with cyclic edge stabilizers. This provides a striking parallel with the decomposition theory of 3-manifolds, a parallel that extends far beyond a mere analogy and accounts for much of the deepest work in geometric group theory in the past ten years. If you want to learn more about this, search the literature for JSJ decompositions. You may also want to search for complexes of groups, which provide the appropriate higher-dimensional analogue for graphs of groups.

# 8 Geometric Group Theory

Let us refresh the image of $K ( P )$ in our mind’s eye by thinking again about the presentation $P \equiv \langle a , b \ |$ $a b a ^ { - 1 } b ^ { - 1 } \rangle$ of $\mathbb { Z } .$ The complex $K ( P )$ , as we saw earlier, is a torus. Now the torus can be defined as the quotient of the Euclidean plane $\mathbb { R } ^ { 2 }$ by the action of the group $\mathbb { Z } ^ { 2 }$ (where the point $( m , n ) \in \mathbb { Z } ^ { 2 }$ acts as the translation $( x , y ) \mapsto ( x + m , y + n ) )$ : in fact, $\mathbb { R } ^ { 2 }$ , with an appropriate square tiling, is the universal cover of the torus. If we look at the orbit of the point 0 under this action, it forms a copy of $ { \mathbb { Z } } ^ { 2 }$ , and one can thereby see the large-scale geometry of $ { \mathbb { Z } } ^ { 2 }$ laid out for us. We can make the idea of the “geometry of $\mathbb { Z } ^ { 2 ^ { , , } }$ precise by decreeing that edges of the tiling have length 1 and defining the graph distance between vertices to be the length of the shortest path of edges connecting them.

As this example shows, the construction of $K ( P )$ involves the two main (intertwined) strands of geometric group theory. In the first and more classical strand, one studies actions of groups on metric and topological spaces in order to elucidate the structures of both the space and the group (as with the action o $\therefore \mathbb { Z } ^ { 2 }$ on the plane in our example, or the action of the fundamental group of $K ( P )$ on its universal cover in general). The quality of the insights that one obtains varies according to whether the action has or does not have certain desirable properties. The action of $ { \mathbb { Z } } ^ { 2 }$ on $\mathbb { R } ^ { 2 }$ consists of isometries on a space with a fine geometric structure, and the quotient (the torus) is compact. Such actions are in many ways ideal, but sometimes one accepts weaker admission criteria in order to obtain a more diverse class of groups, and sometimes one demands even more structure in order to narrow the focus and study groups and spaces of an exceptional, but for that reason interesting, character.

This first strand of geometric group theory mingles with the second. In the second strand, one regards finitely generated groups as geometric objects in their own right equipped with word metrics, which are defined as follows. Given a finite generating set S for a group $T ,$ one defines the Cayley graph of Γ by joining each element $\gamma \in T$ by an edge to each element of the form $\gamma s$ or $\gamma s ^ { - 1 }$ with $s \in S$ (which is the same as the graph formed by the edges of the universal covering of $K ( P ) )$ . The distance $d _ { S } ( \gamma _ { 1 } , \gamma _ { 2 } )$ between $\gamma _ { 1 }$ and $\gamma _ { 2 }$ is then the length of the shortest path from $\gamma _ { 1 }$ to γ2 if all edges have length 1. Equivalently, it is the length of the shortest word in the free group on S that is equal to $\chi _ { 1 } ^ { - 1 } \chi _ { 2 }$ in Γ .

The word metric and the Cayley graph depend on the choice of generating set but their large-scale geometry does not. In order to make this idea precise, we introduce the notion of a quasi-isometry. This is an equivalence relation that identifies spaces that are similar on a large scale. If X and $Y$ are two metric spaces, then a quasi-isometry from X to Y is a function $\phi : X \to Y$ with the following two properties. First, there are positive constants $c , C ,$ and  such that $c d ( x , x ^ { \prime } ) - \epsilon \leqslant$ $d ( \phi ( x ) , \phi ( x ^ { \prime } ) ) \leqslant C d ( x , x ^ { \prime } )$ : this says that $\phi$ distorts sufficiently large distances by at most a constant factor. Second, there is a constant $C ^ { \prime }$ such that for every $y \in Y$ there is some $x \in X$ for which $d ( \phi ( x ) , y ) \leqslant C ^ { \prime } \colon$ : this says that φ is a “quasi-surjection” in the sense that every element of Y is close to the image of an element of X .

Consider for example the two spaces $\mathbb { R } ^ { 2 }$ and $ { \mathbb { Z } } ^ { 2 }$ , where the metric on $ { \mathbb { Z } } ^ { 2 }$ is given by the graph distance defined earlier. In this case the map $\phi : \mathbb { R } ^ { 2 } \ \to \ \mathbb { Z } ^ { 2 }$ that takes $( x , y )$ to $( \lfloor x \rfloor , \lfloor y \rfloor )$ (where x denotes the largest integer less than or equal to x) is easily seen to be a quasi-isometry: if the Euclidean distance d between two points $( x , y )$ and $( x ^ { \prime } , y ^ { \prime } )$ is at least 10, say, then the graph distance between $( \lfloor x \rfloor , \lfloor y \rfloor )$ and $( \lfloor x ^ { \prime } \rfloor , \lfloor y ^ { \prime } \rfloor )$ will certainly lie between ${ \scriptstyle { \frac { 1 } { 2 } } } d$ and 2d. Notice how little we care about the local structure of the two spaces: the map $\phi$ is a quasi-isometry despite not even being continuous.

It is not hard to check that if $\phi$ is a quasi-isometry from X to Y , then there is a quasi-isometry ψ from Y to X that “quasi-inverts” φ, in the sense that every x in X is at most a bounded distance from $\psi \phi ( x )$ and every y in Y is at most a bounded distance from $\phi \psi ( \boldsymbol { y } )$ . Once one has established this, it is easy to see that quasiisometry is an equivalence relation.

Returning to Cayley graphs and word metrics, it turns out that if you take two different sets of generators for the same group, then the resulting Cayley graphs will be quasi-isometric. Thus, any property of a Cayley graph that is invariant under quasi-isometry will be a property not just of the graph but of the group itself. When dealing with such invariants we are free to think of Γ itself as a space (since we do not care which Cayley graph we form), and we can replace it by any metric space that is quasi-isometric to it, such as the universal cover of a closed Riemannian manifold with fundamental group Γ (whose existence we discussed earlier). Then the tools of analysis can be brought to bear on it.

A fundamental fact, discovered independently by many people and often called the Milnor–Švarc lemma, provides a crucial link between the two main strands of geometric group theory. Let us call a metric space X a length space if the distance between each pair of points is the infimum of the lengths of paths joining them. The Milnor–Švarc lemma states that if a group Γ acts “properly discontinuously” as a set of isometries of a length space X, and if the quotient is compact, then Γ is finitely generated and quasi-isometric to X (for any choice of word metric).

We have seen an example of this already: $ { \mathbb { Z } } ^ { 2 }$ is quasiisometric to the Euclidean plane. Less obviously, the same is true of ${ \varGamma } _ { \Delta }$ . (Consider the map that takes each element α of $T _ { \Delta }$ to the point of $ { \mathbb { Z } } ^ { 2 }$ nearest α(0).)

The fundamental group of a compact Riemannian manifold is quasi-isometric to the universal cover of that manifold. Therefore, from the point of view of quasi-isometry invariants, the study of such manifolds is equivalent to the study of arbitrary finitely presented groups. In a moment we will discuss some nontrivial consequences of this equivalence. But first let us reflect on the fact that, when finitely generated groups are considered as metric objects in the framework of large-scale geometry, they present us with a new challenge: we should classify finitely generated groups up to quasi-isometry.

This is an impossible task, of course, but nevertheless serves as a beacon in modern geometric group theory, one that has guided us toward many beautiful theorems, particularly under the general heading of rigidity. For example, suppose that you come across a finitely generated group Γ that is reminiscent of $\mathbb { Z } ^ { n }$ on a large scale: in other words, quasi-isometric to it. We are not necessarily given any algebraically defined map between this mystery group and $\mathbb { Z } ^ { n }$ , and yet it transpires that such a group must contain a copy of $\mathbb { Z } ^ { n }$ as a subgroup of finite index.

At the heart of this result is Gromov’s polynomialgrowth theorem, a landmark theorem published in 1981. This theorem concerns the number of points within a distance r of the identity in a finitely generated group Γ . This will be a function $f ( r )$ , and Gromov was interested in how the function f (r ) grows as r tends to infinity, and what that tells us about the group Γ .

If Γ is an Abelian group with d generators, then it is not hard to see that $f ( r )$ is at most $( 2 r + 1 ) ^ { d }$ (since each generator is raised to a power between r and $r ) .$ Thus, in this case $f ( r )$ is bounded above by a polynomial in r . At the other extreme, if Γ is a free group with two generators a and b, say, then $f ( r )$ is exponentially large, since all sequences of length r that consist of as and bs (and not their inverses) give different elements of Γ .

Given this sharp contrast in behavior, one might wonder whether requiring $f ( r )$ to be bounded above by a polynomial forces Γ to exhibit a great deal of commutativity. Fortunately, there is a much-studied definition that makes this idea precise. Given any group G and any subgroup H of $G ,$ the commutator [G, H] is the subgroup generated by all elements of the form $g h g ^ { - 1 } h ^ { - 1 }$ , where $_ g$ belongs to G and h belongs to H. If G is Abelian, then [G, H] contains just the identity. If G is not Abelian, then [G, G] forms a group $G _ { 1 }$ that contains other elements besides the identity, but it may be that $[ G , G _ { 1 } ]$ is trivial. In that case, one says that G is a two-step nilpotent group. In general, a k-step nilpotent group G is one where, if you form a sequence by setting $G _ { 0 } = G$ and $G _ { i + 1 } = [ G , G _ { i } ]$ for each i, then you eventually reach the trivial group, and the first time you do so is at $G _ { k }$ . A nilpotent group is a group that is k-step nilpotent for some k.

Gromov’s theorem states that a group has polynomial growth if and only if it has a nilpotent subgroup of finite index. This is a quite extraordinary fact: the polynomial-growth condition is easily seen to be independent of the choice of word metric and to be an invariant of quasi-isometry. Thus the seemingly rigid and purely algebraic condition of having a nilpotent subgroup of finite index is in fact a quasi-isometry invariant, and therefore a flabby, robust characteristic of the group.

In the past fifteen years quasi-isometric rigidity theorems have been established for many other classes of groups, including lattices in semisimple Lie groups and the fundamental groups of compact 3-manifolds (where the classification up to quasi-isometry involves more than algebraic equivalences), as well as various classes defined in terms of their graph of group decompositions. In order to prove theorems of this type, one must identify nontrivial invariants of quasi-isometry that allow one to distinguish and relate various classes of spaces. In many cases such invariants come from the development of suitable analogues of the tools of algebraic topology, modified so that they behave well with respect to quasi-isometries rather than continuous maps.

# 9 The Geometry of the Word Problem

It is time to explain the comments I made earlier about the geometry inherent in the basic decision problems of combinatorial group theory. I shall concentrate exclusively on the geometry of the word problem.

Gromov’s filling theorem describes a startlingly intimate connection between the highly geometric study of disks with minimal area in riemannian geometry [I.3 §6.10] and the study of word problems, which seems to belong more to algebra and logic.

On the geometric side, the basic object of study is the isoperimetric function $\mathrm { F i l l } _ { M } ( l )$ of a complete Riemannian manifold M. Given any contractible closed path of length $l ,$ there is a disk of minimal area that is bounded by that path. The largest such area, over all closed paths of length $l ,$ is defined to be $\mathrm { F i l l } _ { M } ( l )$ . Thus, the isoperimetric function is the smallest function of which it is true to say that every closed path of length l can be filled by a disk of area at most $\mathrm { F i l l } _ { M } ( l )$ .

The image to have in mind here is that of a soap film: if one twists a circular wire of length l in Euclidean space and dips it in soap, the film that forms has area at most $l ^ { 2 } / 4 \pi ,$ , whereas if one performs the same experiment in hyperbolic space [I.3 §6.6], the area of the film is bounded by a linear function of l. Correspondingly, the isoperimetric functions of $\mathbb { E } ^ { n }$ and $\mathbb { H } ^ { n }$ (and quotients of them by groups of isometries) are quadratic and linear, respectively. In a moment we shall discuss what types of isoperimetric functions arise when one considers other geometries (more precisely, compact Riemannian manifolds).

To state the filling theorem we need to think about the algebraic side as well. Here, we identify a function that measures the complexity of a direct attack on the word problem for an arbitrary finitely presented group ${ \cal { T } } ~ = ~ \langle { A } ~ | ~ R \rangle$ . If we wish to know whether a word w equals the identity in Γ and do not have any further insight into the nature of Γ , then there is not much we can do other than repeatedly insert or remove the given relations $r \in R$ .

Consider the simple example ${ \cal { T } } = \langle a , b \mid b ^ { 2 } a ,$ , baba . In this group $a b a ^ { 2 } b$ represents the identity. How do we prove this? Well,

$$
\begin{array}{l} a b a ^ {2} b = a (b ^ {2} a) b a ^ {2} b = a b (b a b a) a b \\ = a b a b = a (b a b a) a ^ {- 1} = a a ^ {- 1} = 1. \\ \end{array}
$$

Now let us think about the proof geometrically, via the Cayley graph. Since $a b a ^ { 2 } b ~ = ~ 1$ in the group $T ,$ we obtain a cycle in this graph if we start at the identity and go along edges labeled $a , b , a , a , b ,$ in that order (in which case we visit the vertices $1 , a , a b , a b a , a b a ^ { 2 }$ , $a b a ^ { 2 } b = 1 )$ . The equalities in the proof can be thought of as a way of “contracting” this cycle down to the identity by means of inserting or deleting small loops: for instance, we could insert $b , a , b ,$ a into the list of edge directions, since baba is a relation, or we could delete a trivial loop of the form $a , a ^ { - 1 }$ . This contraction can be given a more topological character if we turn our Cayley graph into a two-dimensional complex by filling in each small loop with a face. Then the contraction of the original cycle consists in gradually moving it across these small faces.

Thus, the difficulty of demonstrating that a word w equals the identity is intimately connected with the area of w, denoted Area(w), which can be thought of algebraically as the smallest sequence of relations you need to insert and delete to turn w into the identity, or geometrically as the smallest number of faces you need to make a disk that fills the cycle represented by w.

The Dehn function $\delta _ { \Gamma } : \mathbb { N } \to \mathbb { N }$ bounds Area(w) in terms of the length w of the word w: $\delta _ { \cal { T } } ( n )$ is the largest area of any word of length at most n that equals 1 in Γ . If the Dehn function grows rapidly, then the word problem is hard, since there are short words that are equal to the identity, but their area is very large, so that any demonstration that they are equal to the identity has to be very long. Results bounding the Dehn function are called isoperimetric inequalities.

The subscript on $\delta _ { \varGamma }$ is somewhat misleading since different finite presentations of the same group will in general yield different Dehn functions. This ambiguity is tolerated because it is tightly controlled: if the groups defined by two finite presentations are isomorphic, or just quasi-isometric, then the corresponding Dehn functions have similar growth rates. More precisely, they are equivalent, with respect to what is sometimes called the standard equivalence relation $" \simeq "$ of geometric group theory: given two monotone functions $f , g : [ 0 , \infty ) \to [ 0 , \infty )$ , one writes $f \preccurlyeq g$ if there exists a constant $C > 0$ such that $f ( l ) \leqslant C g ( C l + C ) + C l + C$ for all $l \geqslant 0$ , and $f \simeq g \operatorname { i f } f \preccurlyeq g$ and $g \preccurlyeq f ;$ and one extends this relation to include functions from N to $[ 0 , \infty )$ .

You will have noticed a resemblance between the definitions of $\mathrm { F i l l } _ { M } ( l )$ and $\delta _ { \cal T } ( n )$ . The filling theorem relates them precisely: it states that if M is a smooth compact manifold, then FillM $( l ) \simeq \delta _ { T } ( l )$ , where Γ is the fundamental group $\pi _ { 1 }$ M of M .

For example, since $ { \mathbb { Z } } ^ { 2 }$ is the fundamental group of the torus $T = \mathbb { R } ^ { 2 } / \mathbb { Z } ^ { 2 }$ , which has Euclidean geometry, $\delta _ { \mathbb { Z } ^ { 2 } } ( l )$ is quadratic.

# 9.1 What Are the Dehn Functions?

We have seen that the complexity of word problems is related to the study of isoperimetric problems in Riemannian and combinatorial geometry. Such insights have, in the last fifteen years, led to great advances in the understanding of the nature of Dehn functions. For example, one can ask for which numbers $\rho$ the function $n ^ { \rho }$ is a Dehn function. The set of all such numbers, which can be shown to be countable, is known as the isoperimetric spectrum, denoted IP, and it is now largely understood.

Following work by many authors, Brady and Bridson proved that the closure of IP is $\{ 1 \} \cup [ 2 , \infty )$ . The finer structure of IP was described by Birget, Rips, and Sapir in terms of the time functions of Turing machines. A further result by the same authors and Ol’shanskii explains how fundamental Dehn functions are to understanding the complexity of arbitrary approaches to the word problem for finitely generated groups Γ : the word problem for Γ lies in NP if and only if Γ is a subgroup of a finitely presented group with polynomial Dehn function. (Here, NP is the class of problems in the famous “P versus $\mathrm { N P ^ { \prime } }$ question: see computational complexity [IV.20 §3] for a description of this class.)

The structure of IP raises an obvious question: What can one say about the two classes of groups singled out as special—those with linear Dehn functions and those with quadratic ones? The true nature of the class of groups with a quadratic Dehn function remains obscure for the moment but there is a beautifully definitive description of those with a linear Dehn function: they are the word hyperbolic groups, which we shall discuss in the next section.

Not all Dehn functions are of the form $n ^ { \alpha } \colon$ there are Dehn functions such as $n ^ { \alpha }$ log n, for example, and others that grow more quickly than any iterated exponential, for example that of

$$
\langle a, b \mid a b a ^ {- 1} b a b ^ {- 1} a ^ {- 1} b ^ {- 2} \rangle .
$$

If the word problem for Γ is unsolvable, then $\delta _ { \cal { T } } ( n )$ will grow faster than any recursive function (indeed this serves as a definition of such groups).

# 9.2 The Word Problem and Geodesics

A closed geodesic on a Riemannian manifold is a loop that locally minimizes distance, such as a loop formed by an elastic band when released on a perfectly smooth surface. Examples such as the great circles on a sphere or the waist of an hourglass show that manifolds may contain closed geodesics that are null-homotopic: that is, they can be moved continuously until they are reduced to a point. But can one construct a compact topological manifold with the property that no matter what metric one puts on it there will always be infinitely many such geodesics? (Technically, if you go around a geodesic loop n times, then you get a geodesic; we avoid this by counting only “primitive” geodesics.)

From a purely geometric point of view this is a daunting problem: all specific metric information has been stripped away and one has to deal with an arbitrary metric on the floppy topological object left behind. But group theory provides a solution: if the Dehn function of the fundamental group π1M grows at least as fast as $2 ^ { 2 ^ { n } }$ , then in any Riemannian metric on M there will be infinitely many closed geodesics that are null-homotopic. The proof of this is too technical to sketch here.

# 10 Which Groups Should One Study?

Several special classes of groups have emerged from our previous discussion, such as nilpotent groups, 3-manifold groups, groups with linear Dehn functions, and groups with a single defining relation. Now we shall change viewpoint and ask which groups present themselves for study as we set out to explore the universe of all finitely presented groups, starting with the easiest ones.

The trivial group comes first, of course, followed by the finite groups. Finite groups are discussed in various other places in this volume, so I shall ignore them in what follows and adopt the approach of large-scale geometry, blurring the distinction between groups that have a common subgroup of finite index.

The first infinite group is surely $\mathbb { Z } ,$ but what comes next is open to debate. If one wants to retain the safety of commutativity, then finitely generated Abelian groups come next. Then, as one slowly relinquishes commutativity and control over growth and constructibility, one passes through the progressively larger classes of nilpotent, polycyclic, solvable, and elementary amenable groups. We have already met nilpotent groups in our discussion of Gromov’s polynomial-growth theorem. They crop up in many contexts as the most natural generalization of Abelian groups and much is known about them, not least because one can prove a great deal by induction on the k for which they are k-step nilpotent. One can also exploit the fact that G is built from the finitely generated Abelian groups $G _ { i } / G _ { i + 1 }$ in a very controlled way. The larger class of polycyclic groups is built in a similar way, while finitely generated solvable groups are built in a finite number of steps from Abelian groups that need not be finitely generated. This last class is not only larger but wilder; the isomorphism problem is solvable among polycyclic groups, for example, but unsolvable among solvable groups. By definition a group G is solvable if its derived series, defined inductively by $G ^ { ( n ) } \ = \ [ G ^ { ( n - 1 ) } , G ^ { ( n - 1 ) } ]$ with $G ^ { ( 0 ) } = G ,$ terminates in a finite number of steps.

The concept known as amenability forms an important link between geometry, analysis, and group theory. Solvable groups are amenable but not vice versa. It is not quite the case that a finitely presented group is amenable if and only if it does not contain a free subgroup of rank 2, but for a novice this serves as a good rule of thumb.

Now, let us return to Z in a more adventurous frame of mind, throw away the security of commutativity, and start taking free products instead. In this more liberated approach, finitely generated free groups appear after Z as the first groups in the universe. What comes next? Thinking geometrically, we might note that free groups are precisely those groups that have a tree as a Cayley graph and then ask which groups have Cayley graphs that are tree-like.

A key property of a tree is that all of its triangles are degenerate: if you take any three points in the tree and join them by shortest paths, then every point in one of these paths is contained in at least one other path as well. This is a manifestation of the fact that trees are spaces of infinite negative curvature. To get a feeling for why, consider what happens when one rescales the metric on a space of bounded negative curvature such as the hyperbolic plane H2. If we replace the standard distance function $d ( x , y )$ by $( 1 / n ) d ( x , y )$ and let n tend $\mathrm { ~ t o ~ } \infty$ , then the curvature of this space (in the classical sense of differential geometry) tends to $- \infty$ . This is captured by the fact that triangles look increasingly degenerate: there is a constant δ(n), with $\delta ( n )  0$ as $n \to \infty$ , such that any side of a triangle in the scaled hyperbolic space $( \mathbb { H } ^ { 2 } , ( 1 / n ) d )$ is contained in the δ(n)-neighborhood of the union of the other two sides. More colloquially, triangles in $\mathbb { H } ^ { 2 }$ are uniformly thin and get increasingly thin as one rescales the metric.

With this picture in mind, one might move a little away from trees by asking which groups have Cayley graphs in which all triangles are uniformly thin. (It makes little sense to specify the thinness constant δ since it will change when one changes generating set.) The answer is Gromov’s hyperbolic groups. This is a fascinating class of groups that has many equivalent definitions and arises in many contexts. For example, we have already met it as the class of groups that have linear Dehn functions. (It is not at all obvious that these two definitions are equivalent.)

Gromov’s great insight is that because the thin-triangles condition encapsulates so much of the essence of the large-scale geometry of negatively curved manifolds, hyperbolic groups share many of the rich properties enjoyed by the groups that act nicely by isometries on such spaces. Thus, for example, hyperbolic groups have only finitely many conjugacy classes of finite subgroups, contain no copy of $\mathbb { Z } ^ { 2 }$ , and (after accounting for torsion) have compact classifying spaces. Their conjugacy problems can be solved in less than quadratic time, and Sela showed that one can even solve the isomorphism problem among torsion-free hyperbolic groups. In addition to their many fascinating properties and natural definition, a further source of interest in hyperbolic groups is the fact that in a precise statistical sense, a random finitely presented group will be hyperbolic.

Spaces of negative and nonpositive curvature have played a central role in many branches of mathematics in the last twenty years. There is no room even to begin to justify this assertion here but it does guide us in where to look for natural enlargements of the class of hyperbolic groups: we want nonpositively curved groups, defined by requiring that their Cayley graphs enjoy a key geometric feature that cocompact groups of isometries inherit from simply connected spaces of nonpositive curvature (“CAT(0) spaces”). But in contrast to the hyperbolic case, the class of groups that one obtains varies considerably when one perturbs the definition, and delineating the resulting classes and their (rich) properties has been the subject of much research.

The added complications that one encounters when one moves from negative to nonpositive curvature are exemplified by the fact that the isomorphism problem is unsolvable in one of the most prominent classes that arises: the so-called combable groups.

Let us now return to free groups and ask which hyperbolic groups are the immediate neighbors of free groups. Remarkably, this vague question has a convincing answer.

One of the great triumphs of arboreal group theory is the proof that there is a finite description of the set Hom(G, F) of homomorphisms from an arbitrary finitely generated group G to a free group F. The basic building blocks in this description are what Sela calls limit groups. One of the many ways of defining a limit group L is that for each finite subset $X \subset L$ there should exist a homomorphism to a finitely generated free group that is injective on X.

Limit groups can also be defined as those whose first-order logic [IV.23 §1] resembles that of a free group in a precise sense. To see how first-order logic can be used to say something nontrivial about a group, consider the sentence

$$
\forall x, y, z
$$

$$
(x y \neq y x) \vee (y z \neq z y) \vee (x z = z x) \vee (y = 1).
$$

A group with this property is commutative transitive: if x commutes with y ≠ 1, and y commutes with z, then x commutes with z. Free groups and Abelian groups have this property but a direct product of non-Abelian free groups, for example, does not.

It is a simple exercise to show that free Abelian groups are limit groups. But if one restricts attention to groups that have precisely the same first-order logic as free groups, one gets a smaller class consisting only of hyperbolic groups. The groups in this class are the subject of intense scrutiny at the moment. They all have negatively curved two-dimensional classifying spaces, built from graphs and hyperbolic surfaces in a hierarchical manner. The fundamental groups $\Sigma _ { g }$ of closed surfaces of genus $g \geqslant 2$ lie in this class, lending substance to the traditional opinion in combinatorial group theory that, among nonfree groups, it is the groups $\Sigma _ { g }$ that resemble free groups $F _ { n }$ most closely.

Incorporating this opinion into our earlier discussion, we arrive at the view that the groups $\mathbb { Z } ^ { n }$ , the free groups $F _ { n } .$ , and the groups $\Sigma _ { g }$ are the most basic of infinite groups. This is the start of a rich vein of ideas involving the automorphisms of these groups. In particular, there are many striking parallels between their outer automorphism groups $\mathrm { G L } _ { n } ( \mathbb { Z } ) , \mathrm { O u t } ( F _ { n } )$ , and $\mathbf { M o d } _ { g } ~ \cong ~ \operatorname { O u t } ( \Sigma _ { g } )$ (the mapping class group). These three classes of groups play a fundamental role across a broad spectrum of mathematics. I have mentioned them here in order to make the point that, beyond the search for knowledge about natural classes of groups, there are certain “gems” in group theory that merit a deep and penetrating study in their own right. Other groups that people might suggest for this category include Coxeter groups (generalized reflection groups, for which ${ \cal { T } } _ { \Delta }$ is a prototype) and Artin groups (particularly braid groups [III.4], which again crop up in many branches of mathematics).

I have thrown classes of groups at you thick and fast in this last section. Even so, there are many fascinating classes of groups and important issues that I have ignored completely. But so it must be, for as Higman’s theorem assures us, the challenges, joys, and frustrations of finitely presented groups can never be exhausted.

# Further Reading

Bridson, M. R., and A. Haefliger. 1999. Metric Spaces of Non-Positive Curvature. Grundlehren der Mathematischen Wissenschaften, volume 319. Berlin: Springer.   
Gromov, M. 1984. Infinite groups as geometric objects. In Proceedings of the International Congress of Mathematicians, Warszawa, Poland, 1983, volume 1, pp. 385–92. Warsaw: PWN.   
. 1993. Asymptotic invariants of infinite groups. In Geometric Group Theory, volume 2. London Mathematical Society Lecture Note Series, volume 182. Cambridge: Cambridge University Press.   
Lyndon, R. C., and P. E. Schupp. 2001. Combinatorial Group Theory. Classics in Mathematics. Berlin: Springer.

# IV.11 Harmonic Analysis

# Terence Tao

# 1 Introduction

Much of analysis tends to revolve around the study of general classes of functions [I.2 §2.2] and operators [III.50]. The functions are often real-valued or complexvalued, but may take values in other sets, such as a vector space [I.3 §2.3] or a manifold [I.3 §6.9]. An operator is itself a function, but at a “second level,” because its domain and range are themselves spaces of functions: that is, an operator takes a function (or perhaps more than one function) as its input and returns a transformed function as its output. Harmonic analysis focuses in particular on the quantitative properties of such functions, and how these quantitative properties change when various operators are applied to them.1

What is a “quantitative property” of a function? Here are two important examples. First, a function is said to be uniformly bounded if there is some real number M such that $| f ( x ) | \leqslant M$ for every x. It can often be useful to know that two functions $f$ and $^ g$ are $^ { \mathrm { 6 } } \mathrm { u n i } \cdot$ formly close,” which means that their difference $f - g$ is uniformly bounded with a small bound M. Second, a function is called square integrable if the integral $\int | f ( x ) | ^ { 2 }$ dx is finite. The square integrable functions are important because they can be analyzed using the theory of hilbert spaces [III.37].

A typical question in harmonic analysis might then be the following: if a function $f : \mathbb { R } ^ { n }$ R is square integrable, its gradient $\nabla f$ exists, and all the n components of $\nabla f$ are also square integrable, does this imply that $f$ is uniformly bounded? (The answer is yes when $n = 1$ , and no, but only just, when $n = 2 ;$ this is a special case of the Sobolev embedding theorem, which is of fundamental importance in the analysis of partial differential equations [IV.12].) If so, what are the precise bounds one can obtain? That is, given the integrals of $| f | ^ { 2 }$ and $| ( \nabla f ) _ { i } | ^ { 2 }$ , what can you say about the uniform bound M that you obtain for $f 2$

Real and complex functions are of course very familiar in mathematics, and one meets them in high school. In many cases one deals primarily with special functions [III.85]: polynomials, exponentials, trigonometric functions, and other very concrete and explicitly defined functions. Such functions typically have a very rich algebraic and geometric structure, and many questions about them can be answered exactly using techniques from algebra and geometry.

However, in many mathematical contexts one has to deal with functions that are not given by an explicit formula. For example, the solutions to ordinary and partial differential equations often cannot be given in an explicit algebraic form (as a composition of familiar functions such as polynomials, exponential functions [III.25], and trigonometric functions [III.92]). In such cases, how does one think about a function? The answer is to focus on its properties and see what can be deduced from them: even if the solution of a differential equation cannot be described by a useful formula, one may well be able to establish certain basic facts about it and be able to derive interesting consequences from those facts. Some examples of properties that one might look at are measurability, boundedness, continu-$\mathrm { i t y } ,$ , differentiability, smoothness, analyticity, integrability, or quick decay at infinity. One is thus led to consider interesting general classes of functions: to form such a class one chooses a property and takes the set of all functions with that property. Generally speaking, analysis is much more concerned with these general classes of functions than with individual functions. (See also function spaces [III.29].)

This approach can in fact be useful even when one is analyzing a single function that is very structured and has an explicit formula. It is not always easy, or even possible, to exploit this structure and formula in a purely algebraic manner, and then one must rely (at least in part) on more analytical tools instead. A typical example is the Airy function

$$
\operatorname{Ai} (x) = \int_ {- \infty} ^ {\infty} \mathrm{e} ^ {\mathrm{i} (x \xi + \xi^ {3})} \mathrm{d} \xi .
$$

Although this is defined explicitly as a certain integral, if one wants to answer such basic questions as whether Ai(x) is always a convergent integral, and whether this integral goes to zero as x  , it is easiest to proceed using the tools of harmonic analysis. In this case, one can use a technique known as the principle of stationary phase to answer both these questions affirmatively, although there is the rather surprising fact that the Airy function decays almost exponentially fast as $x \to + \infty ,$ but only polynomially fast as $x \to - \infty$ .

Harmonic analysis, as a subfield of analysis, is particularly concerned not just with qualitative properties like the ones mentioned earlier, but also with quantitative bounds that relate to those properties. For instance, instead of merely knowing that a function $f$ is bounded, one may wish to know how bounded it is. That is, what is the smallest $M \geqslant 0$ such that $| f ( x ) | \leqslant M$ for all (or almost all) $x \in \mathbb { R } ;$ this number is known as the sup norm or $L ^ { \infty }  – n o r m$ of $f ,$ and is denoted $\| f \| _ { L ^ { \infty } }$ . Or instead of assuming that $f$ is square integrable one can quantify this by introducing the $L ^ { 2 }$ -norm $\begin{array} { r } { \| f \| _ { L ^ { 2 } } = ( \int | f ( x ) | ^ { 2 } \mathrm { d } x ) ^ { 1 / 2 } } \end{array}$ ; more generally one can quantify pth-power integrability for $0 < p < \infty$ via the $L ^ { p } { \cdot } n o r m$ $\begin{array} { r } { \| f \| _ { L ^ { p } } = ( \int | f ( x ) | ^ { p } \mathrm { d } x ) ^ { 1 / p } } \end{array}$ . Similarly, most of the other qualitative properties mentioned above can be quantified by a variety of norms [III.62], which assign a nonnegative number (or ) to any given function and which provide some quantitative measure of one characteristic of that function. Besides being of importance in pure harmonic analysis, quantitative estimates involving these norms are also useful in applied mathematics, for instance in performing an error analysis of some numerical algorithm.

Functions tend to have infinitely many degrees of freedom, and it is thus unsurprising that the number of norms one can place on a function is infinite as well: there are many ways of quantifying how “large” a function is. These norms can often differ quite dramatically from each other. For instance, if a function $f$ is very large for just a few values, so that its graph has tall, thin “spikes,” then it will have a very large $L ^ { \infty }$ -norm, but $\int \mid f ( x ) \mid \mathrm { d } x ,$ its $L ^ { 1 }$ -norm, may well be quite small. Conversely, if f has a very broad and spread-out graph, then it is possible for $\int \mid f ( x )$ dx to be very large even if $\mid f ( x ) \mid$ is small for every x: such a function has a large $L ^ { 1 }$ -norm but a small $L ^ { \infty }$ -norm. Similar examples can be constructed to show that the $L ^ { 2 }$ -norm sometimes behaves very differently from either the $L ^ { 1 }$ -norm or the $L ^ { \infty }$ -norm. However, it turns out that the $L ^ { 2 } .$ -norm lies “between” these two norms, in the sense that if one controls both the $L ^ { 1 }$ -norm and the $L ^ { \infty }$ -norm, then one also automatically controls the $L ^ { 2 }$ -norm. Intuitively, the reason is that if the $L ^ { \infty }$ -norm is not too large then one eliminates all the spiky functions, and if the $L ^ { 1 }$ -norm is small then one eliminates most of the broad functions; the remaining functions end up being well-behaved in the intermediate $L ^ { 2 } .$ -norm. More quantitatively, we have the inequality

$$
\| f \| _ {L ^ {2}} \leqslant \| f \| _ {L ^ {1}} ^ {1 / 2} \| f \| _ {L ^ {\infty}} ^ {1 / 2},
$$

which follows easily from the trivial algebraic fact that if $| f ( x ) | \leqslant M$ , then $| f ( x ) | ^ { 2 } \leqslant M | f ( x )$ |. This inequality is a special case of hölder’s inequality [V.19], which is one of the fundamental inequalities in harmonic analysis. The idea that control of two “extreme” norms automatically implies further control on “intermediate” norms can be generalized tremendously and leads to very powerful and convenient methods known as interpolation, which is another basic tool in this area.

The study of a single function and all its norms eventually gets somewhat tiresome, though. Nearly all fields of mathematics become a lot more interesting when one considers not just objects, but also maps between objects. In our case, the objects in question are functions, and, as was mentioned in the introduction, a map that takes functions to functions is usually referred to as an operator. (In some contexts it is also called a transform [III.91].) Operators may seem like fairly complicated mathematical objects—their inputs and outputs are functions, which in turn have inputs and outputs that are usually numbers—but they are in fact a very natural concept since there are many situations where one wants to transform functions. For example, differentiation can be thought of as an operator, which takes a function $f$ to its derivative ${ \mathrm { d } } f / { \mathrm { d } } x .$ This operator has a well-known (partial) inverse, integration, which takes $f$ to the function F that is defined by the formula

$$
F (x) = \int_ {- \infty} ^ {x} f (y) d y.
$$

A less intuitive, but particularly important, example is the fourier transform [III.27]. This takes $f$ to a function $\hat { f } ,$ given by the formula

$$
\hat {f} (x) = \int_ {- \infty} ^ {\infty} \mathrm{e} ^ {- 2 \pi \mathrm{i} x y} f (y) \mathrm{d} y.
$$

It is also of interest to consider operators that take two or more inputs. Two particularly common examples are the pointwise product and convolution. If $f$ and $_ g$ are two functions, then their pointwise product $f g$ is defined in the obvious way:

$$
(f g) (x) = f (x) g (x).
$$

The convolution, denoted $f * g ,$ is defined as follows:

$$
f * g (x) = \int_ {- \infty} ^ {\infty} f (y) g (x - y) d y.
$$

This is just a very small sample of interesting operators that one might look at. The original purpose of harmonic analysis was to understand the operators that were connected to Fourier analysis, real analysis, and complex analysis. Nowadays, however, the subject has grown considerably, and the methods of harmonic analysis have been brought to bear on a much broader set of operators. For example, they have been particularly fruitful in understanding the solutions of various linear and nonlinear partial differential equations, since the solution of any such equation can be viewed as an operator applied to the initial conditions. They are also very useful in analytic and combinatorial number theory, when one is faced with understanding the oscillation present in various expressions such as exponential sums. Harmonic analysis has also been applied to analyze operators that arise in geometric measure theory, probability theory, ergodic theory, numerical analysis, and differential geometry.

A primary concern of harmonic analysis is to obtain both qualitative and quantitative information about the effects of these operators on generic functions. A typical example of a quantitative estimate is the inequality

$$
\| f * g \| _ {L ^ {\infty}} \leqslant \| f \| _ {L ^ {2}} \| g \| _ {L ^ {2}},
$$

which is true for all $f , g \in L ^ { 2 }$ . This result, which is a special case of Young’s inequality, is easy to prove: one just writes out the definition of $f * g ( x )$ and applies the cauchy–schwarz inequality [V.19]. As a consequence, one can draw the qualitative conclusion that the convolution of two functions in $L ^ { 2 }$ is always continuous. Let us briefly sketch the argument, since it is an instructive one.

A fundamental fact about functions in $L ^ { 2 }$ is that any such function $f$ can be approximated arbitrarily well (in the $L ^ { 2 }$ -norm) by a function $\tilde { f }$ that is continuous and compactly supported. (The second condition means that $\tilde { f }$ takes the value zero everywhere outside some interval [−M, M].) Given any two functions $f$ and $g \mathrm { i n } L ^ { 2 }$ , let $\tilde { f }$ and g˜ be approximations of this kind. It is an exercise in real analysis to prove that $\tilde { f } * \tilde { g }$ is continuous, and it follows easily from the inequality above that $\tilde { f } * \tilde { g }$ is close to $f * g$ in the $L ^ { \infty } .$ -norm, since

$$
f * g - \tilde {f} * \tilde {g} = f * (g - \tilde {g}) + (f - \tilde {f}) * \tilde {g}.
$$

Therefore, $f$ g can be approximated arbitrarily well in the L∞-norm by continuous functions. A standard result in basic real analysis (that a uniform limit of continuous functions is continuous) now tells us that $f * g$ is continuous.

Notice the general structure of this argument, which occurs frequently in harmonic analysis. First, one identifies a “simple” class of functions for which one can easily prove the result one wants. Next, one proves that every function in a much wider class can be approximated in a suitable sense by simple functions. Finally, one uses this information to deduce that the result holds for functions in the wider class as well. In our case, the simple functions were the continuous functions of finite support, the wider class consisted of square-integrable functions, and the suitable sense of approximation was closeness in the $L ^ { 2 } .$ -norm.

We shall give some further examples of qualitative and quantitative analysis of operators in the next section.

# 2 Example: Fourier Summation

To illustrate the interplay between quantitative and qualitative results, we shall now sketch some of the basic theory of summation of Fourier series, which historically was one of the main motivations for studying harmonic analysis.

In this section, we shall consider functions $f$ that are periodic with period 2π: that is, functions such that $f ( x + 2 \pi ) = f ( x )$ for all x. An example of such a function is $f ( x ) = 3 + \sin ( x ) - 2 \cos ( 3 x )$ . A function like this, which can be written as a finite linear combination of functions of the form sin(nx) and cos(nx), is called a trigonometric polynomial. The word “polynomial” is used here because any such function can be expressed as a polynomial in sin(x) and cos(x), or alternatively, and somewhat more conveniently, as a polynomial in $\mathrm { e } ^ { \mathrm { i } x }$ and $\mathbf { e } ^ { - \mathrm { i } x }$ . That is, it can be written as $\textstyle \sum _ { n = - N } ^ { N } c _ { n } \mathrm { e } ^ { \mathrm { i } n x }$ for some N and some choice of coefficients $( c _ { n } : - N \leqslant n \leqslant N )$ . If we know that $f$ can be expressed in this form, then we can work out the coefficient $c _ { n }$ quite easily: it is given by the formula

$$
c _ {n} = \frac {1}{2 \pi} \int_ {0} ^ {2 \pi} f (x) \mathrm{e} ^ {- \mathrm{i} n x} \mathrm{d} x.
$$

It is a remarkable and very important fact that we can say something similar about a much wider class of functions—if, that is, we now allow infinite linear combinations. Suppose that f is a periodic function that is also continuous (or, more generally, that $f$ is absolutely integrable, meaning that the integral of $\mid f ( x ) \mid$ between 0 and $2 \pi$ is finite). We can then define the Fourier coefficients ${ \hat { f } } ( n )$ of $f ,$ , using exactly the formula we had above for $c _ { n } \colon$

$$
\hat {f} (n) = \frac {1}{2 \pi} \int_ {0} ^ {2 \pi} f (x) \mathrm{e} ^ {- \mathrm{i} n x} \mathrm{d} x.
$$

The example of trigonometric polynomials now suggests that one should have the identity

$$
f (x) = \sum_ {n = - \infty} ^ {\infty} \hat {f} (n) \mathrm{e} ^ {\mathrm{i} n x},
$$

expressing $f$ as a sort of “infinite trigonometric polynomial,” but this is not always true, and even when it is true it takes some effort to justify it rigorously, or even to say precisely what the infinite sum means.

To make the question more precise, let us introduce for each natural number N the Dirichlet summation operator $S _ { N }$ . This takes a function f to the function $S _ { N } f$ that is defined by the formula

$$
S _ {N} f (x) = \sum_ {n = - N} ^ {N} \hat {f} (n) \mathrm{e} ^ {\mathrm{i} n x}.
$$

The question we would like to answer is whether $S _ { N } f$ converges to $f$ as $N \ \to \ \infty$ . The answer turns out to be surprisingly complicated: not only does it depend on the assumptions that one places on the function $f ,$ , but it also depends critically on how one defines $^ { * } { \mathrm { c o n } } \cdot$ vergence.” For example, if we assume that $f$ is continuous and ask for the convergence to be uniform, then the answer is very definitely no: there are examples of continuous functions $f$ for which $S _ { N } f$ does not even converge pointwise to $f .$ . However, if we ask for a weaker form of convergence, the answer is yes: $S _ { N } f$ will necessarily converge to $f$ in the $L ^ { p }$ topology for any $0 ~ < ~ p ~ < ~ \infty$ , and even though it does not have to converge pointwise, it will converge almost everywhere, meaning that the set of x for which $S _ { N } f ( x )$ does not converge to x has measure [III.55] zero. If instead one assumes only that $f$ is absolutely integrable, then it is possible for the partial sums $S _ { N } f$ to diverge at every single point $x ,$ as well as being divergent in the $L ^ { p }$ topology for every $p$ such that $0 < p \leqslant \infty$ . The proofs of all of these results ultimately rely on very quantitative results in harmonic analysis, and in particular on various Lp-type estimates on the Dirichlet sum $S _ { N } f ( x )$ , as well as estimates connected with the closely related maximal operator, which takes $f$ to the function $\mathrm { s u p } _ { N > 0 } | S _ { N } f ( x )$ |.

As these results are a little tricky to prove, let us first discuss a simpler result, in which the Dirichlet summation operators $S _ { N }$ are replaced by the Fejér summation operators $F _ { N }$ . For each $N _ { : }$ , the operator $F _ { N }$ is the average of the first N Dirichlet operators: that is, it is given by the formula

$$
F _ {N} = \frac {1}{N} (S _ {0} + \dots + S _ {N - 1}).
$$

It is not hard to show that if $S _ { N } f$ converges to $f ,$ then so does $F _ { N } f$ . However, by averaging the $S _ { N } f$ we allow cancellations to take place that sometimes make it possible for $F _ { N } f$ to converge to $f$ even when $S _ { N } f$ does not. Indeed, here is a sketch of a proof that $F _ { N } f$ converges to f whenever f is continuous and periodic—which, as we have seen, is far from true of $S _ { N } f$ .

In its basic structure, the argument is similar to the one we used when showing that the convolution of two functions in $L ^ { 2 }$ is continuous. Note first that the result is easy to prove when $f$ is a trigonometric polynomial, since then $S _ { N } f \ = \ f$ for every N from some point onward. Now the Weierstrass approximation theorem says that every continuous periodic function $f$ can be uniformly approximated by trigonometric polynomials: that ${ \mathrm { i } } s ,$ for every $\varepsilon > 0$ there is a trigonometric polynomial such that $\| f - g \| _ { L ^ { \infty } } \leqslant \varepsilon .$ We know that $F _ { N } g$ is close to $_ g$ for large N (since $^ g$ is a trigonometric polynomial), and would like to deduce the same for $f .$ The first step is to use some routine trigonometric manipulation to prove the identity

$$
F _ {N} f (x) = \int_ {- \pi} ^ {\pi} \frac {\sin^ {2} (\frac {1}{2} N y)}{N \sin^ {2} (\frac {1}{2} y)} f (x - y) d y.
$$

The precise form of this expression is less important than two properties of the function

$$
u (y) = \frac {\sin^ {2} (\frac {1}{2} N y)}{N \sin^ {2} (\frac {1}{2} y)}
$$

that we shall use. One is that $u ( y )$ is always nonnegative and the other is that $\textstyle \int _ { - \pi } ^ { \pi } u ( y ) \mathrm { d } y = 1$ . These two facts allow us to say that

$$
\begin{array}{l} F _ {N} h (x) = \int_ {- \pi} ^ {\pi} u (y) h (x - y) d y \\ \leqslant \| h \| _ {L ^ {\infty}} \int_ {- \pi} ^ {\pi} u (y) d y = \| h \| _ {L ^ {\infty}}. \\ \end{array}
$$

That is, $\| F _ { N } h \| _ { L ^ { \infty } } \leqslant \| h \| _ { L ^ { \infty } }$ for any bounded function h.

To apply this result, we choose a trigonometric polynomial $_ g$ such that $\| f - g \| _ { L ^ { \infty } } \leqslant \varepsilon$ and let $h = f - g .$ . Then we find that $\| F _ { N } h \| _ { L ^ { \infty } } = \| F _ { N } f - F _ { N } g \| _ { L ^ { \infty } } \leqslant \varepsilon$ as well. As mentioned above, if we choose N large enough, then $\| F _ { N } g - g \| _ { L ^ { \infty } } \leqslant \varepsilon ,$ and then we use the triangle inequality [V.19] to say that

$$
\begin{array}{l} \| F _ {N} f - f \| _ {L ^ {\infty}} \\ \leqslant \left\| F _ {N} f - F _ {N} g \right\| _ {L ^ {\infty}} + \left\| F _ {N} g - g \right\| _ {L ^ {\infty}} + \left\| g - f \right\| _ {L ^ {\infty}}. \\ \end{array}
$$

Since each term on the right-hand side is at most $\varepsilon ,$ this shows that $\| F _ { N } f - f \| _ { L ^ { \infty } }$ is at most 3ε. And since ε can be made arbitrarily small, this shows that $F _ { N } f$ converges to $f .$ .

A similar argument (using minkowski’s integral inequality [V.19] instead of the triangle inequality) shows that $\| F _ { N } f \| _ { L ^ { p } } \leqslant \| f \| _ { L ^ { p } }$ for all $1 ~ \leqslant ~ p ~ \leqslant ~ \infty$ , $f \in L ^ { p }$ , and $N \geqslant 1 .$ . As a consequence, one can modify the above argument to show that $F _ { N } f$ converges to $f$ in the $L ^ { p }$ topology for every $f \in L ^ { p }$ . A slightly more difficult result (relying on a basic result in harmonic analysis known as the Hardy–Littlewood maximal inequality) asserts that, for every $1 < p \leqslant \infty$ , there exists a constant $C _ { p }$ such that one has the inequality $\| \operatorname { s u p } _ { N } | F _ { N } f | \| _ { L ^ { p } } \leqslant C _ { p } \| f \| _ { L ^ { p } }$ for all $f \in L ^ { p } ;$ ; as a consequence, one can show that $F _ { N } f$ converges to $f$ almost everywhere for every $f ~ \in ~ L ^ { p }$ and $1 ~ < ~ p ~ \leqslant ~ \infty .$ A slight modification of this argument also allows one to treat the endpoint case when $f$ is merely assumed to be absolutely integrable; see the discussion on the Hardy–Littlewood maximal inequality at the end of this article.

Now let us return briefly to Dirichlet summation. Using some fairly sophisticated techniques in harmonic analysis (such as Calderón–Zygmund theory) one can show that when $1 < p < \infty$ the Dirichlet operators $S _ { N }$ are bounded in $L ^ { p }$ uniformly in N. In other words, for every $p$ in this range there exists a positive real number $C _ { p }$ such that $\| S _ { N } f \| _ { L ^ { p } } \leqslant C _ { p } \| f \| _ { L ^ { p } }$ for every function $f$ in $L ^ { p }$ and every nonnegative integer $N .$ As a consequence, one can show that $S _ { N } f$ converges to $f$ in the $L ^ { p }$ topology for all $f$ in $L ^ { p }$ and every $p$ such that $1 ~ < ~ p ~ < ~ \infty$ . However, the quantitative estimate on $S _ { N }$ fails at the endpoints $\ p \quad = \ 1$ and $p \quad = \quad \infty ,$ and from this one can also show that the convergence result also fails at these endpoints (either by explicitly constructing a counterexample or by using general results such as the so-called uniform boundedness principle).

What happens if we ask for $S _ { N } f$ to converge to $f$ almost everywhere? Almost-everywhere convergence does not follow from convergence in $L ^ { p }$ when $p \_ <$ < , so we cannot use the above results to prove it. It turns out to be a much harder question, and was a famous open problem, eventually answered by carleson’s theorem [V.5] and an extension of it by Hunt. Carleson proved that one has an estimate of the form sup $\ L _ { \ V } | S _ { N } f | \| _ { L ^ { p } } \leqslant C _ { p } \| f \| _ { L ^ { p } }$ in the case $p \ = \ 2$ , and Hunt generalized the proof to cover all $p$ with $1 \ <$ < $p \ < \ \infty$ . This result implies that the Dirichlet sums of an $L ^ { p }$ function do indeed converge almost everywhere when $1 ~ < ~ p ~ \leqslant ~ \infty$ . On the other hand, this estimate fails at the endpoint $p ~ = ~ 1$ , and there is in fact an example due to kolmogorov [VI.88] of an absolutely integrable function whose Dirichlet sums are everywhere divergent. These results require a lot of harmonic analysis theory. In particular they use many decompositions of both the spatial variable and the frequency variable, keeping the Heisenberg uncertainty principle in mind. They then carefully reassemble the pieces, exploiting various manifestations of orthogonality.

To summarize, quantitative estimates such as $L ^ { p }$ estimates on various operators provide an important route to establishing qualitative results, such as convergence of certain series or sequences. In fact there are a number of principles (notably the uniform boundedness principle and a result known as Stein’s maximal principle) which assert that in certain circumstances this is the only route, in the sense that a quantitative estimate must exist in order for the qualitative result to be true.

# 3 Some General Themes in Harmonic Analysis: Decomposition, Oscillation, and Geometry

One feature of harmonic analysis methods is that they tend to be local rather than global. For instance, if one is analyzing a function f it is quite common to decompose it as a sum $f = f _ { 1 } + \cdot \cdot \cdot + f _ { k } ,$ , with each function $f _ { i }$ “localized” in the sense that its support (the set of values $_ x$ for which $f _ { i } ( x ) \neq 0 )$ has a small diameter. This would be called localization in the spatial variable. One can also localize in the frequency variable by applying the process to the Fourier transform $\hat { \boldsymbol f }$ of $f .$ Having split $f$ up like this, one can carry out estimates for the pieces separately and then recombine them later. One reason for this “divide and conquer” strategy is that a typical function $f$ tends to have many different features—for example, it may be very “spiky,” “discontinuous,” or “high frequency” in some places, and “smooth” or “low frequency” in others—and it is difficult to treat all of these features at once. A well-chosen decomposition of the function f can isolate these features from each other, so that each component has only one salient feature that could cause difficulty: the spiky part can go into one $f _ { i } ,$ the high-frequency part into another, and so on. In reassembling the estimates from the individual components, one can use crude tools such as the triangle inequality or more refined tools, for instance those relying on some sort of orthogonality, or perhaps a clever algorithm that groups the components into manageable clusters. The main drawback of the decomposition method (other than an aesthetic one) is that it tends to give bounds that are not quite optimal; however, in many cases one is content with an estimate that differs from the best possible one by a multiplicative constant.

To give a simple example of the method of decomposition, let us consider the Fourier transform $\hat { f } ( \xi )$ of a function $f : \mathbb { R } \to \mathbb { C } ,$ , defined (for suitably nice functions $f )$ by the formula

$$
\hat {f} (\xi) = \int_ {\mathbb {R}} f (x) \mathrm{e} ^ {- 2 \pi \mathrm{i} x \xi} \mathrm{d} x.
$$

What we can say about the size of $\hat { \boldsymbol f } ,$ as measured by suitable norms, if we are given information about the size of $f ,$ as measured by other norms?

Here are two simple observations in response to this question. First, since the modulus of $\mathbf { e } ^ { - 2 \pi \mathrm { i } x \xi }$ is always equal to 1, it follows that $| { \hat { f } } ( \xi )$ is at most $\int _ { \mathbb { R } } | f ( x ) |$ dx. This tells us that $\| \hat { f } \| _ { L ^ { \infty } } \leqslant \| f \| _ { L ^ { 1 } }$ , at least if $f \in L ^ { 1 }$ . In particular, ${ \hat { f } } \in L ^ { \infty }$ . Secondly, the Plancherel theorem, a very basic fact of Fourier analysis, tells us that $\| \hat { f } \| _ { L ^ { 2 } } \mathrm { i s }$ equal to $\| f \| _ { L ^ { 2 } } { \mathrm { ~ i f ~ } } f \in L ^ { 2 }$ . Therefore, if f belongs to $L ^ { 2 }$ then so does $\hat { f } .$ .

We would now like to know what happens if f lies in an intermediate $L ^ { p }$ space. In other words, what happens if $1 < p < 2 ?$ Since $L ^ { p }$ is not contained in either $L ^ { 1 }$ or $L ^ { 2 }$ , one cannot use either of the above two results directly. However, let us take a function $f \in L ^ { p }$ and consider what the difficulty is. The reason $f$ may not lie in $L ^ { 1 }$ is that it may decay too slowly: for instance, the function $f ( x ) = ( 1 + | x | ) ^ { - 3 / 4 }$ tends to zero more slowly than $1 / x$ as $x \ \to \ \infty$ , so its integral is infinite. However, if we raise $f$ to the power $3 / 2$ we obtain the function $( 1 + | x | ) ^ { - 9 / 8 }$ which decays quickly enough to have a finite integral, so $f$ does belong to $L ^ { 3 / 2 }$ . Similar examples show that the reason $f$ may fail to belong to $L ^ { 2 }$ is that it can have places where it tends to infinity slowly enough for the integral of $| f | ^ { p }$ to be finite but not slowly enough for the integral of $| f | ^ { 2 }$ to be finite.

Notice that these two reasons are completely different. Therefore, we can try to decompose $f$ into two pieces, one consisting of the part where $f$ is large and the other consisting of the part where f is small. That is, we can choose some threshold λ and define $f _ { 1 } ( x )$ to be $f ( x )$ when $| f ( x ) | < \lambda$ and 0 otherwise, and define $f _ { 2 } ( x )$ to be $f ( x )$ when $| f ( x ) | \geqslant \lambda$ and 0 otherwise. Then $f _ { 1 } + f _ { 2 } = f ,$ , and $f _ { 1 }$ and $f _ { 2 }$ are the “small part” and “large part” of $f ,$ respectively.

Because $| f _ { 1 } ( x ) | < \lambda$ for every x, we find that

$$
\left| f _ {1} (x) \right| ^ {2} = \left| f _ {1} (x) \right| ^ {2 - p} \left| f _ {1} (x) \right| ^ {p} <   \lambda^ {2 - p} \left| f _ {1} (x) \right| ^ {p}.
$$

Therefore, $f _ { 1 }$ belongs to $L ^ { 2 }$ and $\| f _ { 1 } \| _ { L ^ { 2 } } \leqslant \lambda ^ { 2 - p } \| f _ { 1 } \| _ { L ^ { p } }$ . Similarly, because $| f _ { 2 } ( x ) | \geqslant \lambda$ whenever $f _ { 2 } ( x ) \neq 0$ , we have the inequality $| f _ { 2 } ( x ) | \leqslant | f _ { 2 } ( x ) | ^ { p } / \lambda ^ { p - 1 }$ for every x, which tells us that $f _ { 2 }$ belongs to $L ^ { 1 }$ and that $\| f _ { 2 } \| _ { L ^ { 1 } } \leqslant$ $\| f _ { 2 } \| _ { L ^ { p } } / \lambda ^ { p - 1 }$ .

From our knowledge about the $L ^ { 2 }$ -norm of $f _ { 1 }$ and the L1-norm of $f _ { 2 }$ we can obtain upper bounds for the $L ^ { 2 }$ -norm of $\ddot { f } _ { 1 }$ and the $L ^ { \infty }$ -norm of ${ \hat { f } } _ { 2 }$ , by our remarks above. By using this strategy for every λ and combining the results in a clever way, one can obtain the Hausdorff–Young inequality, which is the following assertion. Let $p$ lie between 1 and 2 and let $p ^ { \prime }$ be the dual exponent of $\nu ,$ which is the number $p / ( p - 1 )$ . Then there is a constant $C _ { p }$ such that, for every function $f \in L ^ { p }$ , one has the inequality $\| \hat { f } \| _ { L ^ { p ^ { \prime } } } \leqslant C _ { p } \| f \| _ { L ^ { p } } .$ . The particular decomposition method we have used to obtain this result is formally known as the method of real interpolation. It does not give the best possible value of $C _ { p }$ , which turns out to be $p ^ { 1 / 2 p } / ( p ^ { \prime } ) ^ { 1 / 2 p ^ { \prime } }$ , but that requires more delicate methods.

Another basic theme in harmonic analysis is the attempt to quantify the elusive phenomenon of oscillation. Intuitively, if an expression oscillates wildly, then we expect its average value to be relatively small in magnitude, since the positive and negative parts, or in the complex case the parts with a wide range of different arguments, will cancel out. For instance, if a 2π- periodic function $f$ is smooth, then for large n the Fourier coefficient

$$
\hat {f} (n) = \frac {1}{2 \pi} \int_ {- \pi} ^ {\pi} f (x) \mathrm{e} ^ {- \mathrm{i} n x}
$$

will be very small since $\int _ { - \pi } ^ { \pi } \operatorname { e } ^ { - \mathrm { i } n x } = 0$ and the comparatively slow variation in $f ( x )$ is not enough to stop the cancellation occurring. This assertion can easily be proved rigorously by repeated integration by parts. Generalizations of this phenomenon include the so-called principle $o f$ stationary phase, which among other things allows one to obtain precise control on the Airy function Ai(x) discussed earlier. It also yields the Heisenberg uncertainty principle, which relates the decay and smoothness of a function to the decay and smoothness of its Fourier transform.

A somewhat different manifestation of oscillation lies in the principle that if one has a sequence of functions that oscillate in different ways, then their sum should be significantly smaller than the bound that follows from the triangle inequality. Again, this is the result of cancellation that is simply not noticed by the triangle inequality. For instance, the Plancherel theorem in Fourier analysis implies, among other things, that a trigonometric polynomial $\begin{array} { r } { \sum _ { n = - N } ^ { N } c _ { n } \mathrm { e } ^ { \mathrm { i } n x } } \end{array}$ has an $L ^ { 2 }$ -norm of

$$
\left(\frac {1}{2 \pi} \int_ {0} ^ {2 \pi} \left| \sum_ {n = - N} ^ {N} c _ {n} \mathrm{e} ^ {\mathrm{i} n x} \right| ^ {2}\right) ^ {1 / 2} = \left(\sum_ {n = - N} ^ {N} | c _ {n} | ^ {2}\right) ^ {1 / 2}.
$$

This bound (which can also be proved by direct calculation) is smaller than the upper bound of $\textstyle \sum _ { n = - N } ^ { N } | c _ { n } |$ that would be obtained if we simply applied the triangle inequality to the functions $c _ { n } \mathrm { e } ^ { \mathrm { i } n x }$ . This identity can be viewed as a special case of the Pythagorean theorem, together with the observation that the harmonics $\mathrm { e } ^ { \mathrm { i } n x }$ are all orthogonal to each other with respect to the inner product [III.37]

$$
\langle f, g \rangle = \frac {1}{2 \pi} \int_ {0} ^ {2 \pi} f (x) \overline {{{{g (x)}}}} d x.
$$

This concept of orthogonality has been generalized in a number of ways. For instance, there is a more general and robust concept of “almost orthogonality,” which roughly speaking means that the inner products of a collection of functions are small but not necessarily 0.

Many arguments in harmonic analysis will, at some point, involve a combinatorial statement about certain types of geometric objects such as cubes, balls, or boxes. For instance, one useful such statement is the Vitali covering lemma, which asserts that, given any collection $B _ { 1 } , \ldots , B _ { k }$ of balls in Euclidean space $\mathbb { R } ^ { n }$ , there will be a subcollection $B _ { i _ { 1 } } , \ldots , B _ { i _ { m } }$ of balls that are disjoint, but that nevertheless contain a significant fraction of the volume covered by the original balls. To be precise, one can choose the disjoint balls so that

$$
\operatorname{vol} \left(\bigcup_ {j = 1} ^ {m} B _ {i _ {j}}\right) \geqslant 5 ^ {- n} \operatorname{vol} \left(\bigcup_ {j = 1} ^ {k} B _ {j}\right).
$$

(The constant $5 ^ { - n }$ can be improved, but this will not concern us here.) This result is obtained by a “greedy algorithm”: one picks balls one by one, at each stage choosing the largest ball among the $B _ { j }$ that is disjoint from all the balls already selected.

One consequence of the Vitali covering lemma is the Hardy–Littlewood maximal inequality, which we will briefly describe. Given any function $f \in L ^ { 1 } ( \mathbb { R } ^ { n } )$ , any $x \in \mathbb { R } ^ { n }$ , and any $r > 0 ,$ we can calculate the average of $\vert f \vert$ in the n-dimensional sphere $B ( x , r )$ of center x and radius r . Next, we can define the maximal function F of f by letting F(x) be the largest of all these averages as r ranges over all positive real numbers. (More precisely, one takes the supremum.) Then, for each positive real number λ one can define a set $X _ { \lambda }$ to be the set of all x such that $F ( x ) > \lambda$ . The Hardy–Littlewood maximal inequality asserts that the volume of $X _ { \lambda }$ is at most $5 ^ { n } \| f \| _ { L ^ { 1 } } / \lambda . ^ { 2 }$

To prove it, one observes that $X _ { \lambda }$ can be covered by balls B(x, r ) on each of which the integral of f  is at least λ vol(B(x, r )). To this collection of balls one can then apply the Vitali covering lemma, and the result follows. The Hardy–Littlewood maximal inequality is a quantitative result, but it has as a qualitative consequence the Lebesgue differentiation theorem, which asserts the following. If f is any absolutely integrable function defined on Rn, then for almost every $x \in \mathbb { R } ^ { n }$ the averages

$$
\frac {1}{\operatorname{vol} (B (x , r))} \int_ {B (x, r)} f (y) \mathrm{d} y
$$

of $f$ over the Euclidean balls about x tend to $f ( x )$ as $r \quad  \quad 0$ . This example demonstrates the importance of the underlying geometry (in this case, the combinatorics of metric balls) in harmonic analysis.

# Further Reading

Stein, E. M. 1970. Singular Integrals and Differentiability Properties of Functions. Princeton, NJ: Princeton University Press.   
. 1993. Harmonic Analysis. Princeton, NJ: Princeton University Press.   
Wolff, T. H. 2003. Lectures on Harmonic Analysis, edited by I. Łaba and C. Shubin. University Lecture Series, volume 29. Providence, RI: American Mathematical Society.

# IV.12 Partial Differential Equations Sergiu Klainerman

# Introduction

Partial differential equations (or PDEs) are an important class of functional equations: they are equations, or systems of equations, in which the unknowns are functions of more than one variable. As a very crude analogy, PDEs are to functions as polynomial equations (such as $x ^ { 2 } + y ^ { 2 } = 1$ , for example) are to numbers. The distinguishing feature of PDEs, as opposed to more general functional equations, is that they involve not only unknown functions, but also various partial derivatives of those functions, in algebraic combination with each other and with other, fixed, functions. Other important kinds of functional equations are integral equations, which involve various integrals of the unknown functions, and ordinary differential equations (ODEs), in which the unknown functions depend on only one independent variable (such as a time variable t) and the equation involves only ordinary derivatives ${ \mathrm { d } } / { \mathrm { d } } t , { \mathrm { d } } ^ { 2 } / { \mathrm { d } } t ^ { 2 } , { \mathrm { d } } ^ { 3 } / { \mathrm { d } } t ^ { 3 } , \ldots$ . of these functions.

Given the immense scope of the subject the best I can hope to do is to give a very crude perspective on some of the main issues and an even cruder idea of the multitude of current research directions. The difficulty one faces in trying to describe the subject of PDEs starts with its very definition. Is it a unified area of mathematics, devoted to the study of a clearly defined set of objects (in the way that algebraic geometry studies solutions of polynomial equations or topology studies manifolds, for example), or is it rather a collection of separate fields, such as general relativity, several complex variables, or hydrodynamics, each one vast in its own right and centered on a particular, very difficult, equation or class of equations? I will attempt to argue below that, even though there are fundamental difficulties in formulating a general theory of PDEs, one can nevertheless find a remarkable unity between various branches of mathematics and physics that are centered on individual PDEs or classes of PDEs. In particular, certain ideas and methods in PDEs have turned out to be extraordinarily effective across the boundaries of these separate fields. It is thus no surprise that the most successful book ever written about PDEs did not mention PDEs in its title: it was Methods of Mathematical Physics by courant [VI.83] and hilbert [VI.63].

As it is impossible to do full justice to such a huge subject in such limited space I have been forced to leave out many topics and relevant details; in particular, I have said very little about the fundamental issue of breakdown of solutions, and there is no discussion of the main open problems in PDEs. A longer and more detailed version of the article, which includes these topics, can be found at

http://press.princeton.edu/titles/8350.html

# 1 Basic Definitions and Examples

The simplest example of a PDE is the laplace equation [I.3 §5.4]

$$
\Delta u = 0. \tag {1}
$$

Here, Δ is the Laplacian, that is, the differential operator that transforms functions $u \ = \ u ( x _ { 1 } , x _ { 2 } , x _ { 3 } )$ defined from $\mathbb { R } ^ { 3 }$ to R according to the rule

$$
\begin{array}{l} \Delta u (x _ {1}, x _ {2}, x _ {3}) \\ = \partial_ {1} ^ {2} u (x _ {1}, x _ {2}, x _ {3}) + \partial_ {2} ^ {2} u (x _ {1}, x _ {2}, x _ {3}) + \partial_ {3} ^ {2} u (x _ {1}, x _ {2}, x _ {3}), \\ \end{array}
$$

where $\partial _ { 1 } , \partial _ { 2 } , \partial _ { 3 }$ are standard shorthand for the partial derivatives $\partial / \partial x _ { 1 } , \partial / \partial x _ { 2 } , \partial / \partial x _ { 3 }$ . (We will use this shorthand throughout the article.) Two other fundamental examples (also described in [I.3 §5.4]) are the heat equation and the wave equation:

$$
- \partial_ {t} u + k \Delta u = 0, \tag {2}
$$

$$
- \partial_ {t} ^ {2} u + c ^ {2} \Delta u = 0. \tag {3}
$$

In each case one is asked to find a function u that satisfies the corresponding equations. For the Laplace equation u will depend on $x _ { 1 } , x _ { 2 }$ , and $x _ { 3 } ,$ , and for the other two it will depend on t as well. Observe that equations (2) and (3) again involve the symbol $\Delta ,$ but also partial derivatives with respect to the time variable t. The constants k (which is positive) and c are fixed and represent the rate of diffusion and the speed of light, respectively. However, from a mathematical point of view they are not important, since if $u ( t , x _ { 1 } , x _ { 2 } , x _ { 3 } )$ is a solution of (3), for example, then $\nu ( t , x _ { 1 } , x _ { 2 } , x _ { 3 } ) ~ =$ $u ( t , x _ { 1 } / c , x _ { 2 } / c , x _ { 3 } / c )$ satisfies the same equation with $c = 1$ . Thus, when one is studying the equations one can set these constants to be 1. Both equations are called evolution equations because they are supposed to describe the change of a particular physical object as the time parameter t varies. Observe that (1) can be interpreted as a particular case of both (2) and (3): if $u = u ( t , x _ { 1 } , x _ { 2 } , x _ { 3 } )$ is a solution of either (2) or (3) that is independent of t, then $\partial _ { t } u = 0$ , so u must satisfy (1).

In all three examples mentioned above, we tacitly assume that the solutions we are looking for are sufficiently differentiable for the equations to make sense. As we shall see later, one of the important developments in the theory of PDEs was the study of more refined notions of solutions, such as distributions [III.18], which require only weak versions of differentiability.

Here are some further examples of important PDEs. The first is the schrödinger equation [III.83],

$$
\mathrm{i} \partial_ {t} u + k \Delta u = 0, \tag {4}
$$

where u is a function from $\mathbb { R } \times \mathbb { R } ^ { 3 }$ to C. This equation describes the quantum evolution of a massive particle, $k = \hbar / 2 m$ , where $\hbar > 0$ is Planck’s constant and m is the mass of the particle. As with the heat equation, one can set k to equal 1 after a simple change of variables. Though the equation is formally very similar to the heat equation, it has very different qualitative behavior. This illustrates an important general point about PDEs: that small changes in the form of an equation can lead to very different properties of solutions.

A further example is the Klein–Gordon equation

$$
- \partial_ {t} ^ {2} u + c ^ {2} \Delta u - \left(\frac {m c ^ {2}}{\hbar}\right) ^ {2} u = 0. \tag {5}
$$

This is the relativistic counterpart to the Schrödinger equation: the parameter m has the physical interpretation of mass and $m c ^ { 2 }$ has the physical interpretation of rest energy (reflecting Einstein’s famous equation $E = m c ^ { 2 } )$ . One can normalize the constants c and $m c ^ { 2 } / \hbar$ so that they both equal 1 by applying a suitable change of variables to time and space.

Though all five equations mentioned above first appeared in connection with specific physical phenomena, such as heat transfer for (2) and propagation of electromagnetic waves for (3), they have, miraculously, a range of relevance far beyond their original applications. In particular there is no reason to restrict their study to three space dimensions: it is very easy to generalize them to similar equations in n variables $x _ { 1 } , x _ { 2 } , \ldots , x _ { n }$ .

All the PDEs listed so far obey a simple but fundamental property called the principle of superposition: if $u _ { 1 }$ and $u _ { 2 }$ are two solutions to one of these equations, then any linear combination $\alpha _ { 1 } u _ { 1 }$ +a2u2 of these solutions is also a solution. In other words, the space of all solutions is a vector space [I.3 §2.3]. Equations that obey this property are known as homogeneous linear equations. If the space of solutions is an affine space (that is, a translate of a vector space) rather than a vector space, we say that the PDE is an inhomogeneous linear equation; a good example is Poisson’s equation:

$$
\Delta u = f, \tag {6}
$$

where $f : \mathbb { R } ^ { 3 }  \mathbb { R }$ is a function that is given to us and $u : \mathbb { R } ^ { 3 } $ R is the unknown function. Equations that are neither homogeneous linear nor inhomogeneous linear are known as nonlinear. The following equation, the minimal surface equation [III.94 §3.1], is manifestly nonlinear:

$$
\begin{array}{l} \partial_ {1} \left(\frac {\partial_ {1} u}{(1 + | \partial_ {1} u | ^ {2} + | \partial_ {2} u | ^ {2}) ^ {1 / 2}}\right) \\ + \partial_ {2} \left(\frac {\partial_ {2} u}{(1 + | \partial_ {1} u | ^ {2} + | \partial_ {y} u | ^ {2}) ^ {1 / 2}}\right) = 0. \tag {7} \\ \end{array}
$$

The graphs of solutions $u : \mathbb { R } ^ { 2 } $ R of this equation are area-minimizing surfaces (like soap films).

Equations (1), (2), (3), (4), (5) are not just linear: they are all examples of constant-coefficient linear equations. This means that they can be expressed in the form

$$
\mathcal {P} [ u ] = 0, \tag {8}
$$

where is a differential operator that involves linear combinations, with constant real or complex coefficients, of mixed partial derivatives of $u .$ (Such operators are called constant-coefficient linear differential operators.) For instance, in the case of the Laplace equation (1),  is simply the Laplacian Δ, while for the wave equation (3),  is the d’Alembertian

$$
\mathcal {P} = \Box = - \partial_ {t} ^ {2} + \partial_ {1} ^ {2} + \partial_ {2} ^ {2} + \partial_ {3} ^ {2}.
$$

The characteristic feature of linear constant-coefficient operators is translation invariance. Roughly speaking, this means that if you translate a function $u ,$ then you translate $\mathcal { P } u$ in the same way. More precisely, if $\nu ( x )$ is defined to be u(x a) (so the value of u at x becomes the value of v at $x + a ;$ note that x and a belong to $\mathbb { R } ^ { 3 }$ here), then $\mathcal P \nu ( x )$ is equal to $\mathcal { P } u ( x - a )$ . As a consequence of this basic fact we infer that solutions to the homogeneous, linear, constant-coefficient equation (8) are still solutions when translated.

Since symmetries play such a fundamental role in PDEs we should stop for a moment to make a general definition. A symmetry of a PDE is any invertible operation $T : u \mapsto T ( u )$ from functions to functions that preserves the space of solutions, in the sense that u solves the PDE if and only if T (u) solves the same PDE. A PDE with this property is then said to be invariant under the symmetry T . The symmetry T is often a linear operation, though this does not have to be the case. The composition of two symmetries is again a symmetry, as is the inverse of a symmetry, and so it is natural to view a collection of symmetries as forming a group [I.3 §2.1] (which is typically a finite- or infinite-dimensional lie group [III.48 §1]).

Because the translation group is intimately connected with the fourier transform [III.27] (indeed, the latter can be viewed as the representation theory of the former), this symmetry strongly suggests that Fourier analysis should be a useful tool to solve constant-coefficient PDEs, and this is indeed the case.

Our basic constant-coefficient linear operators, the Laplacian Δ and the d’Alembertian , are formally similar in many respects. The Laplacian is fundamentally associated with the geometry of euclidean space $[ \mathrm { I } . 3 \ S 6 . 2 ] \ \mathbb { R } ^ { 3 }$ and the d’Alembertian is similarly associated with the geometry of minkowski space [I.3 §6.8] $\mathbb { R } ^ { 1 + 3 }$ . This means that the Laplacian commutes with all the rigid motions of the Euclidean space $\mathbb { R } ^ { 3 }$ , while the d’Alembertian commutes with the corresponding class of Poincaré transformations of Minkowski spacetime. In the former case this simply means that invariance applies to all transformations of $\mathbb { R } ^ { 3 }$ that preserve the Euclidean distances between points. In the case of the wave equation, the Euclidean distance has to be replaced by the spacetime distance between points (which would be called events in the language of relativity): if $P = ( t , x _ { 1 } , x _ { 2 } , x _ { 3 } )$ and $Q ( s , y _ { 1 } , y _ { 2 } , y _ { 3 } )$ ), then the distance between them is given by the formula

$$
\begin{array}{l} d _ {M} (P, Q) ^ {2} \\ = - (t - s) ^ {2} + (x _ {1} - y _ {1}) ^ {2} + (x _ {2} - y _ {2}) ^ {2} + (x _ {3} - y _ {3}) ^ {2}. \\ \end{array}
$$

As a consequence of this basic fact we infer that all solutions to the wave equation (3) are invariant under translations and lorentz transformations [I.3 §6.8].

Our other evolution equations (2) and (4) are clearly invariant under rotations of the space variables $x =$ $( x ^ { 1 } , x ^ { 2 } , x ^ { 3 } ) ~ \in ~ \mathbb { R } ^ { 3 }$ , when t is fixed. They are also Galilean invariant, which means, in the particular case of the Schrödinger equation $\left( 4 \right) ,$ , that whenever $u \ =$

1 $\mathbf { \xi } _ { t } ( t , x )$ is a solution so is the function $u _ { \nu } ( t , x ) =$ $\mathrm { e } ^ { \mathrm { i } ( \boldsymbol { x } \cdot \nu ) } \mathrm { e } ^ { \mathrm { i } t | \nu | ^ { 2 } } ( t , \boldsymbol { x } - \nu t )$ for any vector $\nu \in \mathbb { R } ^ { 3 }$ .

Poisson’s equation (6), on the other hand, is an example of a constant-coefficient inhomogeneous linear equation, which means that it takes the form

$$
\mathcal {P} [ u ] = f \tag {9}
$$

for some constant-coefficient linear differential operator and known function f . To solve such an equation requires one to understand the invertibility or otherwise of the linear operator : if it is invertible then u will equal ${ \mathcal { P } } ^ { - 1 } f _ { : }$ , and if it is not invertible then either there will be no solution or there will be infinitely many solutions. Inhomogeneous equations are closely related to their homogeneous counterpart; for instance, if $u _ { 1 } ,$ , $u _ { 2 }$ both solve the inhomogeneous equation (9) with the same inhomogeneous term $f ,$ then their difference $u _ { 1 } - u _ { 2 }$ solves the corresponding homogeneous equation (8).

Linear homogeneous PDEs satisfy the principle of superposition but they do not have to be translation invariant. For example, suppose that we modify the heat equation (2) so that the coefficient k is no longer constant but rather an arbitrary, positive, smooth function of $( x _ { 1 } , x _ { 2 } , x _ { 3 } )$ . Such an equation models the flow of heat in a medium in which the rate of diffusion varies from point to point. The corresponding space of solutions is not translation invariant (which is not surprising as the medium in which the heat flows is not translation invariant). Equations like this are called linear equations with variable coefficients. It is more difficult to solve them and describe their qualitative features than it is for constant-coefficient equations. (See, for example, stochastic processes [IV.24 §5.2] for an approach to equations of type (2) with variable k.) Finally, nonlinear equations such as (7) can often still be written in the form $( 8 ) ,$ but the operator  is now a nonlinear differential operator. For instance, the relevant operator for (7) is given by the formula

$$
\mathcal {P} [ u ] = \sum_ {i = 1} ^ {2} \partial_ {i} \bigg (\frac {1}{(1 + | \partial u | ^ {2}) ^ {1 / 2}} \partial_ {i} u \bigg),
$$

where $| \partial u | ^ { 2 } = ( \partial _ { 1 } u ) ^ { 2 } + ( \partial _ { 2 } u ) ^ { 2 }$ . Operators such as these are clearly not linear. However, because they are ultimately constructed from algebraic operations and partial derivatives, both of which are “local” operations, we observe the important fact that  is at least still a “local” operator. More precisely, if $u _ { 1 }$ and $u _ { 2 }$ are two functions that agree on some open set $D ,$ then the expressions $\mathcal { P } [ u _ { 1 } ]$ and $\mathcal { P } [ u _ { 2 } ]$ also agree on this set. In particular, if $\mathcal { P } [ 0 ] = 0$ (as is the case in our example), then whenever u vanishes on a domain, $\mathcal { P } [ u ]$ will also vanish on that domain.

So far we have tacitly assumed that our equations take place in the whole of a space such as $\mathbb { R } ^ { 3 } , \mathbb { R } ^ { + } \times \mathbb { R } ^ { 3 }$ , or $\mathbb { R } \times \mathbb { R } ^ { 3 }$ . In reality one is often restricted to a fixed domain of that space. Thus, for example, equation (1) is usually studied on a bounded open domain of $\mathbb { R } ^ { 3 }$ subject to a specified boundary condition. Here are some basic examples of boundary conditions.

Example. The Dirichlet problem for Laplace’s equation on an open domain of $D \subset \mathbb { R } ^ { 3 }$ is the problem of finding a function u that behaves in a prescribed way on the boundary of D and obeys the Laplace equation inside.

More precisely, one specifies a continuous function $u _ { 0 } : \partial D \ .$ → R and looks for a continuous function u, defined on the closure $\bar { D }$ of D, that is twice continuously differentiable inside D and solves the equations

$$
\left. \begin{array}{l l} \Delta u (x) = 0 & \text { for   all } x \in D, \\ u (x) = u _ {0} (x) & \text { for   all } x \in \partial D. \end{array} \right\} \tag {10}
$$

A basic result in PDEs asserts that if the domain D has a sufficiently smooth boundary, then there is exactly one solution to the problem (10) for any prescribed function $u _ { 0 }$ on the boundary ∂D.

Example. The Plateau problem is the problem of finding the surface of minimal total area that bounds a given curve.

When the surface is the graph of a function u on some suitably smooth domain $D ,$ in other words a set of the form $\{ ( x , y , u ( x , y ) ) : ( x , y ) \in D \}$ , and the bounding curve is the graph of a function u0 over the boundary ∂D of D, then this problem turns out to be equivalent to the Dirichlet problem (10), but with the linear equation (1) replaced by the nonlinear equation (7). For the above equations, it is also often natural to replace the Dirichlet boundary condition $u ( x ) = u _ { 0 } ( x )$ on the boundary ∂D with another boundary condition, such as the Neumann boundary condition $n ( x ) \cdot \nabla _ { x } u ( x ) = u _ { 1 } ( x )$ on ∂D, where n(x) is the outward normal (of unit length) to D at x. Generally speaking, Dirichlet boundary conditions correspond to “absorbing” or “fixed” barriers in physics, whereas Neumann boundary conditions correspond to “reflecting” or “free” barriers.

Natural boundary conditions can also be imposed for our evolution equations (2)–(4). The simplest one is to prescribe the values of u when t 0. We can think of this more geometrically. We are prescribing the values of u at each spacetime point of form $( 0 , x , y , z )$ , and the set of all such points is a hyperplane in $\mathbb { R } ^ { 1 + 3 }$ : it is an example of an initial time surface.

Example. The Cauchy problem (or initial value problem, sometimes abbreviated to IVP) for the heat equation (2) asks for a solution u $: \mathbb { R } ^ { + } \times \mathbb { R } ^ { 3 } \to \mathbb { R }$ on the spacetime domain $\mathbb { R } ^ { + } \times \mathbb { R } ^ { 3 } = \{ ( t , x ) : t > 0 , \ x \in \mathbb { R } ^ { 3 } \}$ , which equals a prescribed function $u _ { 0 } : \mathbb { R } ^ { 3 }  \mathbb { R }$ on the initial time surface $\{ 0 \} \times \mathbb { R } ^ { 3 } = \partial ( \mathbb { R } ^ { + } \times \mathbb { R } ^ { 3 } )$ .

In other words, the Cauchy problem asks for a sufficiently smooth function $u ,$ defined on the closure of $\mathbb { R } ^ { + } \times \mathbb { R } ^ { 3 }$ and taking values in R, that satisfies the conditions

$$
\left. \begin{array}{l} - \partial_ {t} u (t, x) + k \Delta u (t, x) = 0 \\ \text {   for   every   } (t, x) \in \mathbb {R} ^ {+} \times \mathbb {R} ^ {3}, \\ u (0, x) = u _ {0} (x) \quad \text {   for   every   } x \in \mathbb {R} ^ {3}. \end{array} \right\} \tag {11}
$$

The function $u _ { 0 }$ is often referred to as the initial conditions, or initial data, or just data, for the problem. Under suitable smoothness and decay conditions, one can show that this equation has exactly one solution u for each choice of data $u _ { 0 }$ . Interestingly, this assertion fails if one replaces the future domain $\mathbb { R } ^ { + } \times \mathbb { R } ^ { 3 } =$ $\{ ( t , x ) : t > 0 , x \in \mathbb { R } ^ { 3 } \}$ by the past domain $\mathbb { R } ^ { - } \times \mathbb { R } ^ { 3 } =$ $\{ ( t , x ) : t < 0 , x \in \mathbb { R } ^ { 3 } \}$ .

A similar formulation of the IVP holds for the Schrö- dinger equation (4), though in this case we can solve both to the past and to the future. However, in the case of the wave equation (3) we need to specify not just the initial position $u ( 0 , x ) = u _ { 0 } ( x )$ on the initial time surface $t = 0$ , but also an initial velocity $\partial _ { t } u ( 0 , x ) = u _ { 1 } ( x )$ , since equation (3) (unlike (2) or (4)) cannot formally determine $\partial _ { t } u$ in terms of u. One can construct unique smooth solutions (both to the future and to the past of the initial hyperplane $t = 0 )$ to the IVP for (3) for very general smooth initial conditions $u _ { 0 } , u _ { 1 }$ .

Many other boundary-value problems are possible. For instance, when analyzing the evolution of a wave in a bounded domain D (such as a sound wave), it is natural to work with the spacetime domain $\mathbb { R } \times D$ and prescribe both Cauchy data (on the initial boundary $0 \times D )$ and Dirichlet or Neumann data (on the spatial boundary $\mathbb { R } \times \partial D )$ . On the other hand, when the physical problem under consideration is the evolution of a wave outside a bounded obstacle (for example, an electromagnetic wave), one considers instead the evolution in $\mathbb { R } \times ( \mathbb { R } ^ { 3 } \setminus D )$ with a boundary condition on D.

The choice of boundary condition and initial conditions for a given PDE is very important. For equations of physical interest these arise naturally from the context in which they are derived. For example, in the case of a vibrating string, which is described by solutions of the one-dimensional wave equation $\partial _ { t } ^ { 2 } u - \partial _ { x } ^ { 2 } u = 0$ i n the domain $( a , b ) \times \mathbb { R }$ , the initial conditions $u = u _ { 0 }$ and $\partial _ { t } u = u _ { 1 } \mathrm { a t } t = t _ { 0 }$ amount to specifying the original position and velocity of the string. The boundary condition $u ( a ) = u ( b ) = 0$ is what tells us that the two ends of the string are fixed.

So far we have considered just scalar equations. These are equations where there is only one unknown function u, which takes values either in the real numbers R or in the complex numbers C. However, many important PDEs involve either multiple unknown scalar functions or (equivalently) functions that take values in a multidimensional vector space such as $\mathbb { R } ^ { m }$ . In such cases, we say that we have a system of PDEs. An important example of a system is that of the cauchy– riemann equations [I.3 §5.6]:

$$
\partial_ {1} u _ {2} - \partial_ {2} u _ {1} = 0, \quad \partial_ {1} u _ {1} + \partial_ {2} u _ {2} = 0, \tag {12}
$$

where $u _ { 1 } , u _ { 2 } : \mathbb { R } ^ { 2 }  \mathbb { R }$ are real-valued functions on the plane. It was observed by cauchy [VI.29] that a complex function $w ( x + \mathrm { i } y ) = u _ { 1 } ( x , y ) + \mathrm { i } u _ { 2 } ( x , y )$ is holomorphic $\left[ \mathrm { I . 3 5 5 . 6 } \right]$ if and only if its real and imaginary parts $u _ { 1 } , u _ { 2 }$ satisfy the system (12). This system can still be represented in the form of a constant-coefficient linear PDE (8), but u is now a vector $\textstyle { \bigl ( } { } { } _ { u _ { 2 } } ^ { u _ { 1 } } { \bigr ) }$ , and  is not a scalar differential operator, but rather a matrix of operators $\big ( \begin{array} { c } { { - \hat { \sigma } _ { 2 } \hat { \sigma } _ { 1 } } } \\ { { \hat { \sigma } _ { 1 } \hat { \sigma } _ { 2 } } } \end{array} \big )$ ).

The system (12) contains two equations and two unknowns. This is the standard situation for a determined system. Roughly speaking, a system is called overdetermined if it contains more equations than unknowns and underdetermined if it contains fewer equations than unknowns. Underdetermined equations typically have infinitely many solutions for any given set of prescribed data; conversely, overdetermined equations tend to have no solutions at all, unless some additional compatibility conditions are imposed on the prescribed data.

Observe also that the Cauchy–Riemann operator P has the following remarkable property:

$$
\mathcal {P} ^ {2} [ u ] = \mathcal {P} [ \mathcal {P} [ u ] ] = \binom{\Delta u _ {1}}{\Delta u _ {2}}.
$$

Thus $\mathcal { P }$ can be viewed as a square root of the twodimensional Laplacian Δ. One can define a similar type of square root for the Laplacian in higher dimensions and, more surprisingly, even for the d’Alembertian operator  in $\mathbb { R } ^ { 1 + 3 }$ . To achieve this we need to have four 4 × 4 complex matrices $\boldsymbol { \gamma } ^ { 1 } , \boldsymbol { \gamma } ^ { 2 } , \boldsymbol { \gamma } ^ { 3 } , \boldsymbol { \gamma } ^ { 4 }$ that satisfy the property

$$
\gamma^ {\alpha} \gamma^ {\beta} + \gamma^ {\beta} \gamma^ {\alpha} = - 2 m ^ {\alpha \beta} I.
$$

Here, I is the unit 4  4 matrix and $\begin{array} { r } { m ^ { \alpha \beta } = \frac { 1 } { 2 } } \end{array}$ when $\alpha =$ $\beta = 1 , - \frac { 1 } { 2 }$ when $\alpha = \beta \neq 1$ , and 0 otherwise. Using the γ matrices we can introduce the Dirac operator as follows. If $u ~ = ~ ( u _ { 1 } , u _ { 2 } , u _ { 3 } , u _ { 4 } )$ is a function in $\mathbb { R } ^ { 1 + 3 }$ with values in $\mathbb { C } ^ { 4 }$ , then we set $\mathcal { D } u = \mathrm { i } \gamma ^ { \alpha } \partial _ { \alpha } u .$ . It is easy to check that, indeed, $\mathcal { D } ^ { 2 } u = \ d \Omega u$ . The equation

$$
\mathcal {D} u = k u \tag {13}
$$

is called the Dirac equation and it is associated with a free, massive, relativistic particle such as an electron.

One can extend the concept of a PDE further to cover unknowns that are not, strictly speaking, functions taking values in a vector space, but are instead sections of a vector bundle $[ \mathrm { I V } . 6 \ S 5 ] ,$ or perhaps a map from one manifold [I.3 §6.9] to another; such generalized PDEs play an important role in geometry and modern physics. A fundamental example is given by the einstein field equations [IV.13]. In the simplest, “vacuum,” case, they take the form

$$
\operatorname{Ric} (g) = 0, \tag {14}
$$

where $\operatorname { R i c } ( g )$ is the ricci curvature [III.78] tensor of the spacetime manifold $M \ = \ ( M , g )$ . In this case the spacetime metric itself is the unknown to be solved for. One can often reduce such equations locally to more traditional PDE systems by selecting a suitable choice of coordinates, but the task of selecting a “good” choice of coordinates, and working out how different choices are compatible with each other, is a nontrivial and important one. Indeed, the task of selecting a good set of coordinates in order to solve a PDE can end up being a significant PDE problem in its own right.

PDEs are ubiquitous throughout mathematics and science. They provide the basic mathematical framework for some of the most important physical theories: elasticity, hydrodynamics, electromagnetism, general relativity, and nonrelativistic quantum mechanics, for example. The more modern relativistic quantum field theories lead, in principle, to equations in an infinite number of unknowns, which lie beyond the scope of PDEs. Yet, even in that case, the basic equations preserve the locality property of PDEs. Moreover, the starting point of a quantum field theory [IV.17 §2.1.4] is always a classical field theory, which is described by systems of PDEs. This is the case, for example, in the standard model of weak and strong interactions, which is based on the so-called Yang–Mills–Higgs field theory. If we also include the ordinary differential equations of classical mechanics, which can be viewed as onedimensional PDEs, we see that essentially all of physics is described by differential equations. As examples of PDEs underlying some of our most basic physical theories we refer to the articles that discuss the euler and navier–stokes equations [III.23], the heat equation [III.36], the schrödinger equation [III.83], and the einstein equations [IV.13].

An important feature of the main PDEs is their apparent universality. Thus, for example, the wave equation, first introduced by d’alembert [VI.20] to describe the motion of a vibrating string, was later found to be connected with the propagation of sound and electromagnetic waves. The heat equation, first introduced by fourier [VI.25] to describe heat propagation, appears in many other situations in which dissipative effects play an important role. The same can be said about the Laplace equation, the Schrödinger equation, and many other basic equations.

It is even more surprising that equations that were originally introduced to describe specific physical phenomena have played a fundamental role in several areas of mathematics that are considered to be “pure,” such as complex analysis, differential geometry, topology, and algebraic geometry. Complex analysis, for example, which studies the properties of holomorphic functions, can be regarded as the study of solutions to the Cauchy–Riemann equations (12) in a domain of $\mathbb { R } ^ { 2 }$ . Hodge theory is based on studying the space of solutions to a class of linear systems of PDEs on manifolds that generalize the Cauchy–Riemann equations: it plays a fundamental role in topology and algebraic geometry. the atiyah–singer index theorem [V.2] is formulated in terms of a special class of linear PDEs on manifolds, related to the Euclidean version of the Dirac operator. Important geometric problems can be reduced to finding solutions to specific PDEs, typically nonlinear. We have already seen one example: the Plateau problem of finding surfaces of minimal total area that pass through a given curve. Another striking example is the uniformization theorem [V.34] in the theory of surfaces, which takes a compact Riemannian surface S (a two-dimensional surface with a riemannian metric [I.3 §6.10]) and, by solving the PDE

$$
\Delta_ {S} u + \mathrm{e} ^ {2 u} = K \tag {15}
$$

(which is a nonlinear variant of the Laplace equation (1)), uniformizes the metric so that it is “equally curved” at all points on the surface (or, more precisely, has constant scalar curvature [III.78]) without changing the conformal class of the metric (i.e., without distorting any of the angles subtended by curves on the surface). This theorem is of fundamental importance to the theory of such surfaces: in particular, it allows one to give a topological classification of compact surfaces in terms of a single number χ(S), which is called the euler characteristic [I.4 §2.2] of the surface S. The three-dimensional analogue of the uniformization theorem, the geometrization conjecture [IV.7 §2.4] of Thurston, has recently been established by Perelman, who did so by solving yet another PDE; in this case, the equation is the ricci flow [III.78] equation

$$
\partial_ {t} g = 2 \operatorname{Ric} (g), \tag {16}
$$

which can be transformed into a nonlinear version of the heat equation (2) after a carefully chosen change of coordinates. The proof of the geometrization conjecture is a decisive step toward the total classification of all three-dimensional compact manifolds, in particular establishing the well-known poincaré conjecture [IV.7 §2.4]. To overcome the many technical details in establishing this conjecture, one needs to make a detailed qualitative analysis of the behavior of solutions to the Ricci flow equation, a task which requires just about all the advances made in geometric PDEs in the last hundred years.

Finally, we note that PDEs arise not only in physics and geometry but also in many fields of applied science. In engineering, for example, one often wants to control some feature of the solution u to a PDE by carefully selecting whatever components of the given data one can directly influence; consider, for instance, how a violinist controls the solution to the vibrating string equation (closely related to (3)) by modulating the force and motion of a bow on that string in order to produce a beautiful sound. The mathematical theory dealing with these types of issues is called control theory.

When dealing with complex physical systems, one cannot possibly have complete information about the state of the system at any given time. Instead, one often makes certain randomness assumptions about various factors that influence it. This leads to the very important class of equations called stochastic differential equations (SDEs), where one or more components of the equation involve a random variable [III.71 §4] of some sort. An example of this is in the black–scholes model [VII.9 §2] in mathematical finance. A general discussion of SDEs can be found in stochastic processes [IV.24 §6].

The plan for the rest of this article is as follows. In section 2 I shall describe some of the basic notions and achievements of the general theory of PDEs. The main point I want to make here is that, in contrast with ordinary differential equations, for which a general theory is both possible and useful, partial differential equations do not lend themselves to a useful general theoretical treatment because of some important obstructions that I shall try to describe. One is thus forced to discuss special classes of equations such as elliptic, parabolic, hyperbolic, and dispersive equations. In section 3 I will try to argue that, despite the impossibility of developing a useful general theory that encompasses all, or most, of the important examples, there is nevertheless an impressive unifying body of concepts and methods for dealing with various basic equations, and this gives PDEs the feel of a well-defined area of mathematics. In section 4 I develop this further by trying to identify some common features in the derivation of the main equations that are dealt with in the subject. An additional source of unity for PDEs is the central role played by the issues of regularity and breakdown of solutions, which is discussed only briefly here. In the final section we shall discuss some of the main goals that can be identified as driving the subject.

# 2 General Equations

One might expect, after looking at other areas of mathematics such as algebraic geometry or topology, that there was a very general theory of PDEs that could be specialized to various specific cases. As I shall argue below, this point of view is seriously flawed and very much out of fashion. It does, however, have important merits, which I hope to illustrate in this section. I shall avoid giving formal definitions and focus instead on representative examples. The reader who wants more precise definitions can consult the online version of this article.

For simplicity we shall look mostly at determined systems of PDEs. The simplest distinction, which we have already made, is between scalar equations, such as (1)–(5), which consist of only one equation and one unknown, and systems of equations, such as (12) and (13). Another simple but important concept is that of the order of a PDE, which is defined to be the highest derivative that appears in the equation; this concept is analogous to that of the degree of a polynomial. For instance, the five basic equations (1)–(5) listed earlier are second order in space, although some (such as (2) or (4)) are only first order in time. Equations (12) and (13), as well as the Maxwell equations, are first order.1

We have seen that PDEs can be divided into linear and nonlinear equations, with the linear equations being divided further into constant-coefficient and variablecoefficient equations. One can also divide nonlinear PDEs into several further classes depending on the “strength” of the nonlinearity. At one end of the scale, a semilinear equation is one in which all the nonlinear components of the equation have strictly lower order than the linear components. For instance, equation (15) is semilinear, because the nonlinear component $\mathrm { e } ^ { u }$ is of zero order, i.e., it contains no derivatives, whereas the linear component ΔSu is of second order. These equations are close enough to being linear that they can often be effectively viewed as perturbations of a linear equation. A more strongly nonlinear class of equations is that of quasilinear equations, in which the highestorder derivatives of u appear in the equation only in a linear manner but the coefficients attached to those derivatives may depend in some nonlinear manner on lower-order derivatives. For instance, the second-order equation (7) is quasilinear, because if one uses the product rule to expand the equation, then it takes the quasilinear form

$$
\begin{array}{l} F _ {1 1} \left(\partial_ {1} u, \partial_ {2} u\right) \partial_ {1} ^ {2} u + F _ {1 2} \left(\partial_ {1} u, \partial_ {2} u\right) \partial_ {1} \partial_ {2} u \\ + F _ {2 2} (\partial_ {1} u, \partial_ {2} u) \partial_ {2} ^ {2} u = 0 \\ \end{array}
$$

for some explicit algebraic functions $F _ { 1 1 } , F _ { 1 2 } , F _ { 2 2 }$ of the lower-order derivatives of u. While quasilinear equations can still sometimes be analyzed by perturbative techniques, this is generally more difficult to accomplish than it is for an analogous semilinear equation. Finally, we have fully nonlinear equations, which exhibit no linearity properties whatsoever. A typical example is the Monge–Ampère equation

$$
\det (\mathrm{D} ^ {2} u) = F (x, u, \mathrm{D} u),
$$

where $u : \mathbb { R } ^ { n }  \mathbb { R }$ is the unknown function, Du is the gradient [I.3 §5.3] of $u , \mathrm { D } ^ { 2 } u \mathrm { ~ = ~ } ( \partial _ { i } \partial _ { j } u ) _ { 1 \leqslant i , j \leqslant n }$ is the Hessian matrix of u, and ${ \cal F } : \mathbb { R } ^ { n } \times \mathbb { R } \times \mathbb { R } ^ { n } $ R is a given function. This equation arises in many geometric contexts, ranging from manifold-embedding problems to the complex geometry of calabi–yau manifolds [III.6]. Fully nonlinear equations are among the most difficult and least well-understood of all PDEs.

Remark. Most of the basic equations of physics, such as the Einstein equations, are quasilinear. However, fully nonlinear equations arise in the theory of characteristics of linear PDEs, which we discuss below, and also in geometry.

# 2.1 First-Order Scalar Equations

It turns out that first-order scalar PDEs in any number of dimensions can be reduced to systems of firstorder ODEs. As a simple illustration of this important fact consider the following equation in two space dimensions:

$$
\begin{array}{l} a ^ {1} (x ^ {1}, x ^ {2}) \partial_ {1} u (x ^ {1}, x ^ {2}) + a ^ {2} (x ^ {1}, x ^ {2}) \partial_ {2} u (x ^ {1}, x ^ {2}) \\ = f (x ^ {1}, x ^ {2}), \quad (1 7) \\ \end{array}
$$

where $a ^ { 1 } , a ^ { 2 } , f$ are given real functions in the variables $x = ( x ^ { 1 } , x ^ { 2 } ) \in \mathbb { R } ^ { 2 }$ . We associate with (17) the first-order $2 \times 2$ system

$$
\left. \begin{array}{r l} \frac {\mathrm{d} x ^ {1}}{\mathrm{d} s} (s) & = a ^ {1} (x ^ {1} (s), x ^ {2} (s)), \\ \frac {\mathrm{d} x ^ {2}}{\mathrm{d} s} & = a ^ {2} (x ^ {1} (s), x ^ {2} (s)). \end{array} \right\} \tag {18}
$$

To simplify matters, let us assume that $f = 0 .$ .

Suppose now that $x ( s ) = ( x ^ { 1 } ( s ) , x ^ { 2 } ( s ) )$ is a solution of (18), and let us consider how $u ( x ^ { 1 } ( s ) , x ^ { 2 } ( s ) )$ varies as s varies. By the chain rule we know that

$$
\frac {\mathrm{d}}{\mathrm{d} s} u = \partial_ {1} u \frac {\mathrm{d}}{\mathrm{d} s} \frac {\mathrm{d} x ^ {1}}{\mathrm{d} s} + \partial_ {2} u \frac {\mathrm{d} x ^ {2}}{\mathrm{d} s},
$$

and equations (17) and (18) imply that this equals zero (by our assumption that $f \ = \ 0 ) .$ . In other words, any solution $u \ = \ u ( x ^ { 1 } , x ^ { 2 } )$ of (17) with $f ~ = ~ 0$ is constant along any parametrized curve of the form $x ( s ) =$ $( x ^ { 1 } ( s ) , x ^ { 2 } ( s ) )$ that satisfies (18).

Thus, in principle, if we know the solutions to (18), which are called characteristic curves for the equation (17), then we can find all solutions to (17). I say “in principle” because, in general, the nonlinear system (18) is not so easy to solve. Nevertheless, ODEs are simpler to deal with, and the fundamental theorem of ODEs, which we will discuss later in this section, allows us to solve (18) at least locally and for a small interval in s.

The fact that u is constant along characteristic curves allows us to obtain important qualitative information even when we cannot find explicit solutions. For example, suppose that the coefficients $a ^ { 1 } , a ^ { 2 }$ are smooth (or real analytic) and that the initial data is smooth (or real analytic) everywhere on the set  where it is defined, except at some point $x _ { 0 }$ where it is discontinuous. Then the solution u remains smooth (or real analytic) at all points except along the characteristic curve Γ that starts at x0, or, in other words, along the solution to (18) that satisfies the initial condition $x ( 0 ) = x _ { 0 }$ . That is, the discontinuity at $x _ { 0 }$ propagates precisely along Γ . We see here the simplest manifestation of an important principle, which we shall explain in more detail later: singularities of solutions to PDEs propagate along characteristics (or, more generally, hypersurfaces).

One can generalize equation (17) to allow the coefficients $a _ { 1 } , a _ { 2 } ,$ , and $f$ to depend not only on $x = ( x ^ { 1 } , x ^ { 2 } )$ ) but also on u:

$$
a ^ {1} (x, u (x)) \partial_ {1} u (x) + a ^ {2} (x, u (x)) \partial_ {2} u (x) = f (x, u (x)). \tag {19}
$$

The associated characteristic system becomes

$$
\left. \begin{array}{l} \frac {\mathrm{d} x ^ {1}}{\mathrm{d} s} (s) = a ^ {1} (x (s), u (s, x (s))), \\ \frac {\mathrm{d} x ^ {2}}{\mathrm{d} s} (s) = a ^ {2} (x (s), u (s, x (s))). \end{array} \right\} \tag {20}
$$

As a special example of (19) consider the scalar equation in two space dimensions,

$$
\partial_ {t} u + u \partial_ {x} u = 0, \quad u (0, x) = u _ {0} (x), \tag {21}
$$

which is called the Burgers equation. Here we have set $a ^ { 1 } ( x , u ( x ) ) \ = \ 1$ and $a ^ { 2 } ( x , u ( x ) ) \ = \ u ( x )$ . With this choice of $a ^ { 1 } , a ^ { 2 }$ , we can take $x ^ { 1 } ( s )$ to be s in (20). Then, renaming $x ^ { 2 } ( s )$ as $x ( s )$ , we derive the characteristic equation in the form

$$
\frac {\mathrm{d} x}{\mathrm{d} s} (s) = u (s, x (s)). \tag {22}
$$

For any given solution u of (21) and any characteristic curve $( s , x ( s ) )$ we have $( \mathrm { d } / \mathrm { d } s ) u ( s , x ( s ) ) = 0$ . Thus, in principle, knowing the solutions to (22) should allow us to determine the solutions to (21). However, this argument seems worryingly circular, since u itself appears in (22).

To see how this difficulty can be circumvented, consider the IVP for (21): that is, look for solutions that satisfy $u ( 0 , x ) = u _ { 0 } ( x )$ . Consider an associated characteristic curve x(s) such that, initially, $x ( 0 ) ~ = ~ x _ { 0 }$ . Then, since u is constant along the curve, we must have $u ( s , x ( s ) ) = u _ { 0 } ( x _ { 0 } )$ . Hence, going back to (22), we infer that dx $/ \mathrm { d } s = u _ { 0 } ( x _ { 0 } )$ ) and thus $x ( s ) = x _ { 0 } + s u _ { 0 } ( x _ { 0 } )$ . We thus deduce that

$$
u (s, x _ {0} + s u _ {0} (x _ {0})) = u _ {0} (x _ {0}), \tag {23}
$$

which implicitly gives us the form of the solution u. We see once more, from (23), that if the initial data is smooth (or real analytic) everywhere except at a point $x _ { 0 }$ of the line $t = 0 ,$ , then the corresponding solution is also smooth (or real analytic) everywhere in a small neighborhood V of $x _ { 0 }$ , except along the characteristic curve that begins at $x _ { 0 }$ . The smallness of V is necessary here because new singularities can form at large scales. Indeed, u has to be constant along the lines $x + s u _ { 0 } ( x )$ , whose slopes depend on $u _ { 0 } ( x )$ . At a point where these lines cross we would obtain different values of $u ,$ which is impossible unless u becomes singular by this point. This blow-up phenomenon occurs for any smooth, nonconstant initial data $u _ { 0 }$ .

Remark. There is an important difference between the linear equation (17) and the quasilinear equation (19). The characteristics of the first depend only on the coefficients $a ^ { 1 } ( x ) , a ^ { 2 } ( x )$ , while the characteristics of the second depend explicitly on a particular solution u of the equation. In both cases, singularities can only propagate along the characteristic curves of the equation. For nonlinear equations, however, new singularities can form at large distance scales, whatever the smoothness of the initial data.

The above procedure extends to fully nonlinear scalar equations in $\mathbb { R } ^ { d }$ such as the Hamilton–Jacobi equation

$$
\partial_ {t} u + H (x, \mathrm{D} u) = 0, \quad u (0, x) = u _ {0} (x), \tag {24}
$$

where $u : \mathbb { R } { \times } \mathbb { R } ^ { n } \to \mathbb { R }$ is the unknown function, Du is the gradient of u, and the hamiltonian [III.35] $H : \mathbb { R } ^ { d } \times$ $\mathbb { R } ^ { d }$ R and the initial data $u _ { 0 } : \mathbb { R } ^ { d }$ R are given. For instance, the eikonal equation $\partial _ { t } u \ : = \ : | \mathrm Ḋ u Ḍ Ḍ |$ is a special instance of a Hamilton–Jacobi equation. We associate with (24) the ODE system

$$
\left. \begin{array}{l} \frac {\mathrm{d} x ^ {i}}{\mathrm{d} t} = \frac {\partial}{\partial p _ {i}} H (x (t), p (t)), \\ \frac {\mathrm{d} p _ {i}}{\mathrm{d} t} = - \frac {\partial}{\partial x ^ {i}} H (x (t), p (t)), \end{array} \right\} \tag {25}
$$

where i runs from 1 to d. The equations (25) are known as a Hamiltonian system of ODEs. The relationship between this system and the corresponding Hamilton–Jacobi equation is a little more involved than in the cases discussed above. Briefly, we can construct a solution u to (24) based only on the knowledge of the solutions $( x ( t ) , p ( t ) )$ to (25), which are called the bicharacteristic curves of the nonlinear PDE. Once again, singularities can only propagate along bicharacteristic curves (or hypersurfaces). As in the case of the Burgers equation, singularities will occur for more or less any smooth data. Thus, a classical, continuously differentiable solution can only be constructed locally in time. Both Hamilton–Jacobi equations and Hamiltonian systems play a fundamental role in classical mechanics as well as in the theory of the propagation of singularities in linear PDEs. The deep connection between Hamiltonian systems and first-order Hamilton–Jacobi equations played an important role in the introduction of the Schrödinger equation into quantum mechanics.

# 2.2 The Initial Value Problem for ODEs

Before we can continue with our general presentation of PDEs we need first to discuss, for the sake of comparison, the IVP for ODEs. Let us start with a first-order ODE

$$
\partial_ {x} u (x) = f (x, u (x)) \tag {26}
$$

subject to the initial condition

$$
u (x _ {0}) = u _ {0}. \tag {27}
$$

Let us also assume for simplicity that (26) is a scalar equation and that f is a well-behaved function of x and u, such as $f ( x , u ) = u ^ { 3 } - u + 1$ sin x. From the initial data $u _ { 0 }$ we can determine $\partial _ { x } u ( x _ { 0 } )$ by substituting $x _ { 0 }$ into (26). If we now differentiate the equation (26) with respect to x and apply the chain rule, we derive the equation

$$
\partial_ {x} ^ {2} u (x) = \partial_ {x} f (x, u (x)) + \partial_ {u} f (x, u (x)) \partial_ {x} u (x),
$$

which for the example just defined works out to be $\cos x + 3 u ^ { 2 } ( x ) \partial _ { x } u ( x ) - \partial _ { x } u ( x )$ . Hence,

$$
\partial_ {x} ^ {2} u (x _ {0}) = \partial_ {x} f (x _ {0}, u _ {0}) + \partial_ {u} f (x _ {0}, u _ {0}) \partial_ {x} u _ {0},
$$

and since $\partial _ { x } u ( x _ { 0 } )$ has already been determined we find that $\partial _ { x } ^ { 2 } u ( x _ { 0 } )$ can also be explicitly calculated from the initial data u0. This calculation also involves the function f and its first partial derivatives. Taking higher derivatives of the equation (26) we can recursively determine $\partial _ { x } ^ { 3 } u ( x _ { 0 } )$ , as well as all other higher derivatives of u at x0. Therefore, one can in principle determine u(x) with the help of the Taylor series

$$
\begin{array}{l} u (x) = \sum_ {k \geqslant 0} \frac {1}{k !} \partial_ {x} ^ {k} u \left(x _ {0}\right) \left(x - x _ {0}\right) ^ {k} \\ = u (x _ {0}) + \partial_ {x} u (x _ {0}) (x - x _ {0}) \\ + \frac {1}{2 !} \partial_ {x} ^ {2} (x _ {0}) (x - x _ {0}) ^ {2} + \dots . \\ \end{array}
$$

We say “in principle” because there is no guarantee that the series converges. There is, however, a very important theorem, called the Cauchy–Kovalevskaya theorem, which asserts that if the function f is real analytic, as is certainly the case for our function $f ( x , u ) =$ $u ^ { 3 } - u + 1 +$ sin x, then there will be some neighborhood J of $x _ { 0 }$ where the Taylor series converges to a real-analytic solution u of the equation. It is then easy to show that the solution thus obtained is the unique solution to (26) that satisfies the initial condition (27). To summarize: if f is a well-behaved function, then the initial value problem for ODEs has a solution, at least in some time interval, and that solution is unique.

The same result does not always hold if we consider a more general equation of the form

$$
a (x, u (x)) \partial_ {x} u = f (x, u (x)), \quad u (x _ {0}) = u _ {0}. \tag {28}
$$

Indeed, the recursive argument outlined above breaks down in the case of the scalar equation $( x - x _ { 0 } ) \partial _ { x } u =$ $f ( x , u )$ for the simple reason that we cannot even determine $\partial _ { x } u ( x _ { 0 } )$ from the initial condition $u ( x _ { 0 } ) =$ u0. A similar problem occurs for the equation (u $u _ { 0 } ) \partial _ { x } u = f ( x , u )$ . An obvious condition that allows us to extend our previous recursive argument to (28) is to insist that $a ( x _ { 0 } , u _ { 0 } ) \neq 0$ . Otherwise, we say that the IVP (28) is characteristic. If both a and f are also real analytic, the Cauchy–Kovalevskaya theorem applies again and we obtain a unique, real-analytic solution of (28) in a small neighborhood of $x _ { 0 }$ . In the case of an $N \times N$ system,

$$
A (x, u (x)) \partial_ {x} u = F (x, u (x)), \quad u (x _ {0}) = u _ {0},
$$

$A = A ( x , u )$ is an $N \times N$ matrix, and the noncharacteristic condition becomes

$$
\det A (x _ {0}, u _ {0}) \neq 0. \tag {29}
$$

It turns out, and this is extremely important in the development of the theory of ODEs, that, while the nondegeneracy condition (29) is essential to obtain a unique solution of the equation, the analyticity condition is not at all important: it can be replaced by a simple local Lipschitz condition for A and F. It suffices to assume, for example, that their first partial derivatives exist and that they are locally bounded. This is always the case if the first derivatives of A and F are continuous.

Theorem (the fundamental theorem of ODEs). If the matrix $A ( x _ { 0 } , u _ { 0 } )$ is invertible and if A and F are continuous and have locally bounded first derivatives, then there is some time interval J  R that contains $x _ { 0 } ,$ and a unique solution2 u defined on J that satisfies the initial conditions $u ( x _ { 0 } ) = u _ { 0 }$ .

The proof of the theorem is based on the Picard iteration method. The idea is to construct a sequence of approximate solutions $u _ { ( n ) } ( x )$ that converge to the desired solution. Without loss of generality we can assume A to be the identity matrix.3 One starts by setting $u _ { ( 0 ) } ( x ) = u _ { 0 }$ and then defines, recursively,

$$
\partial_ {x} u _ {(n)} (x) = F (x, u _ {(n - 1)} (x)), \quad u _ {(n - 1)} (x _ {0}) = u _ {0}.
$$

Observe that at every stage all we need to solve is a very simple linear problem, which makes Picard iteration easy to implement numerically. As we shall see below, variations of this method are also used for solving nonlinear PDEs.

Remark. In general, the local existence theorem is sharp, in the sense that its conditions cannot be relaxed. We have seen that the invertibility condition for $A ( x _ { 0 } , u _ { 0 } )$ is necessary. Also, it is not always possible to extend the interval J in which the solution exists to the whole of the real line. As an example, consider the nonlinear equation $\partial _ { x } u \ = \ u ^ { 2 }$ with initial data $u ~ = ~ u _ { 0 }$ at $x \ = \ 0 ,$ for which the solution $u = u _ { 0 } / ( 1 - x u _ { 0 } )$ becomes infinite in finite time: in the terminology of PDEs, it blows up.

In view of the fundamental theorem and the example mentioned above, one can define the main goals of the mathematical theory of ODEs as follows.

(i) Find criteria for global existence. In the case of blow-up describe the limiting behavior.   
(ii) In the case of global existence describe the asymptotic behavior of solutions and families of solutions.

Though it is impossible to develop a general theory that achieves both goals (in practice one is forced to restrict oneself to special classes of equations motivated by applications), the general local existence and uniqueness theorem mentioned above provides a powerful unifying theme. It would be very helpful if a similar situation were to hold for general PDEs.

# 2.3 The Initial Value Problem for PDEs

In the one-dimensional situation one specifies initial conditions at a point. The natural higher-dimensional analogue is to specify them on hypersurfaces $\mathcal { H } \subset \mathbb { R } ^ { d }$ , that is, (d − 1)-dimensional subsets (or, to be more precise, submanifolds). For a general equation of order $k ,$ that is, one that involves k derivatives, we need to specify the values of u and of its first k 1 derivatives in the direction normal to . For example, in the case of the second-order wave equation (3) and the initial hyperplane $t = 0$ we need to specify initial data for u and $\partial _ { t } u .$ .

If we wish to use initial data of this kind to start obtaining a solution, it is important that the data should not be degenerate. (We have already seen this in the case of ODEs.) For this reason, we make the following general definition.

Definition. Suppose that we have a kth-order quasilinear system of equations, and the initial data comes in the form of the first k  1 normal derivatives that a solution u must satisfy on a hypersurface ${ \mathcal { H } } .$ We say that the system is noncharacteristic at a point $x _ { 0 }$ of if we can use the initial data to determine formally all the other higher partial derivatives of u at $x _ { 0 } ,$ , in terms of the data.

As a very rough picture to have in mind, it may be helpful to imagine an infinitesimally small neighborhood of $x _ { 0 } .$ If the hypersurface  is smooth, then its intersection with this neighborhood will be a piece of a (d 1)-dimensional affine subspace. The values of u and the first k 1 normal derivatives on this intersection are given by the initial data, and the problem of determining the other partial derivatives is a problem in linear algebra (because everything is infinitesimally small). To say that the system is noncharacteristic at $x _ { 0 }$ is to say that this linear algebra problem can be uniquely solved, which is the case provided that a certain matrix is invertible. This is the nondegeneracy condition referred to earlier.

To illustrate the idea, let us look at first-order equations in two space dimensions. In this case is a curve $T ,$ and since $k - 1 = 0$ we must specify the restriction of u to $T \subset \mathbb { R } ^ { 2 }$ but we do not have to worry about any derivatives. Thus, we are trying to solve the system

$$
\begin{array}{l} a ^ {1} (x, u (x)) \partial_ {1} u (x) + a ^ {2} (x, u (x)) \partial_ {2} u (x) \\ = f (x, u (x)), \quad u | _ {\Gamma} = u _ {0}, \tag {30} \\ \end{array}
$$

where $a ^ { 1 } , \ a ^ { 2 } ,$ , and $f$ are real-valued functions of $_ x$ (which belongs to $\mathbb { R } ^ { 2 } )$ and u. Assume that in a small neighborhood of a point $p$ the curve Γ is described parametrically as the set of points $\boldsymbol { x } = ( x ^ { 1 } ( s ) , x ^ { 2 } ( s ) )$ . We denote by $n ( s ) = ( n _ { 1 } ( s ) , n _ { 2 } ( s ) )$ a unit normal to Γ .

As in the case of ODEs, which we looked at earlier, we would like to find conditions on Γ such that for a given point in Γ we can determine all derivatives of u from the data $u _ { 0 } ,$ the derivatives of u along Γ , and the equation (30). Out of all possible curves Γ we distinguish in particular the characteristic ones we have already encountered above (see (20)):

$$
\left. \begin{array}{l} \frac {\mathrm{d} x ^ {1}}{\mathrm{d} s} = a ^ {1} (x (s), u (x (s))), \\ \frac {\mathrm{d} x ^ {2}}{\mathrm{d} s} = a ^ {2} (x (s), u (x (s))), \end{array} \right\} \quad x (0) = p.
$$

One can prove the following fact:

Along a characteristic curve, the equation (30) is degenerate. That is, we cannot determine the first-order derivatives of u uniquely in terms of the data u0.

In terms of the rough picture above, at each point there is a direction such that if the hypersurface, which in this case is a line, is along that direction, then the resulting matrix is singular. If you follow this direction, then you travel along a characteristic curve.

Conversely, if the nondegeneracy condition

$$
a ^ {1} (p, u (p)) n _ {1} (p) + a _ {2} (p, u (p)) n _ {2} (p) \neq 0 \tag {31}
$$

is satisfied at some point $p = x ( 0 ) \in T$ , then we can determine all higher derivatives of u at x0 uniquely in terms of the data $u _ { 0 }$ and its derivatives along Γ . If the curve Γ is given by the equation $\psi ( x ^ { 1 } , x ^ { 2 } ) = 0 $ , with nonvanishing gradient $\mathrm { D } \psi ( p ) \neq 0 ,$ then the condition (31) takes the form

$$
a ^ {1} (p, u (p)) \partial_ {1} \psi (p) + a ^ {2} (p, u (p)) \partial_ {2} \psi (p) \neq 0.
$$

With a little more work one can extend the above discussion to higher-order equations in higher dimensions, and even to systems of equations. Particularly important is the case of a second-order scalar equation in $\mathbb { R } ^ { d }$ ,

$$
\sum_ {i, j = 1} ^ {d} a ^ {i j} (x) \partial_ {i} \partial_ {j} u = f (x, u (x)), \tag {32}
$$

together with a hypersurface H in Rd defined by the equation ψ(x)  0, where ψ is a function with nonvanishing gradient Dψ. Define the unit normal at a point $x _ { 0 } ~ \in ~ { \mathcal { H } }$ to be $n \ = \ \mathrm { D } \psi / | \mathrm { D } \psi | ,$ or, in component form, $n _ { i } ~ = ~ \partial _ { i } \psi / | \partial \psi |$ . As initial conditions for (32) we prescribe the values of u and its normal derivative $n [ u ] ( x ) = n _ { 1 } ( x ) \partial _ { 1 } u ( x ) + n _ { 2 } ( x ) \partial _ { 2 } u ( x ) + \cdot \cdot \cdot +$ $n _ { d } ( x ) \partial _ { d } u ( x )$ on :

$$
u (x) = u _ {0} (x), \quad n [ u ] (x) = u _ {1} (x), \quad x \in \mathcal {H}.
$$

It can be shown that  is noncharacteristic (with respect to equation (32)) at a point p (that is, we can determine all derivatives of u at p in terms of the initial data $u _ { 0 } , u _ { 1 } )$ if and only if

$$
\sum_ {i, j = 1} ^ {d} a ^ {i j} (p) \partial_ {i} \psi (p) \partial_ {j} \psi (p) \neq 0. \tag {33}
$$

On the other hand,  is a characteristic hypersurface for (32) if

$$
\sum_ {i, j = 1} ^ {d} a ^ {i j} (x) \partial_ {i} \psi (x) \partial_ {j} \psi (x) = 0 \tag {34}
$$

for every x in H .

Example. If the coefficients a of (32) satisfy the condition

$$
\sum_ {i, j = 1} ^ {d} a ^ {i j} (x) \xi_ {i} \xi_ {j} > 0, \quad \forall \xi \in \mathbb {R} ^ {d}, \forall x \in \mathbb {R} ^ {d}, \tag {35}
$$

then clearly, by (34), no surface in $\mathbb { R } ^ { d }$ can be characteristic. This is the case, in particular, for the Laplace equation Δu  f . Consider also the minimal surface equation (7) written in the form

$$
\sum_ {i, j = 1, 2} h ^ {i j} (\partial u) \partial_ {i} \partial_ {j} u = 0, \tag {36}
$$

with $h ^ { 1 1 } ( \partial u ) \ = \ 1 \ + \ ( \partial _ { 2 } u ) ^ { 2 } , \ h ^ { 2 2 } ( \partial u ) \ = \ 1 \ + \ ( \partial _ { 1 } u ) ^ { 2 } ,$ $h ^ { 1 2 } ( \partial u ) \ = \ h ^ { 2 1 } ( \partial u ) \ = \ - \partial _ { 1 } u \partial _ { 2 } u$ . It is easy to check that the quadratic form associated with the symmetric matrix hij (∂u) is positive definite for every ∂u. Indeed,

$$
h ^ {i j} (\partial u) \xi_ {i} \xi_ {j}
$$

$$
= (1 + | \partial u | ^ {2}) ^ {- 1 / 2} (| \xi | ^ {2} - (1 + | \partial u | ^ {2}) ^ {- 1} (\xi \cdot \partial u) ^ {2}) > 0.
$$

Thus, even though (36) is not linear, we see that all surfaces in $\mathbb { R } ^ { 2 }$ are noncharacteristic.

Example. Consider the wave equation u $\begin{array} { r } { = f \operatorname { i n } \mathbb { R } ^ { 1 + d } . } \end{array}$ All hypersurfaces of the form ψ(t, x) 0 for which

$$
(\partial_ {t} \psi) ^ {2} = \sum_ {i = 1} ^ {d} (\partial_ {i} \psi) ^ {2} \tag {37}
$$

are characteristic. This is the famous eikonal equation, which plays a fundamental role in the study of wave propagation. Observe that it splits into two Hamilton– Jacobi equations (see (24)):

$$
\partial_ {t} \psi = \pm \left(\sum_ {i = 1} ^ {d} (\partial_ {i} \psi) ^ {2}\right) ^ {1 / 2}. \tag {38}
$$

The bicharacteristic curves of the associated Hamiltonians are called bicharacteristic curves of the wave equation. As particular solutions of (37) we find $\psi _ { + } ( t , x ) =$ $( t - t _ { 0 } ) + | x - x _ { 0 } |$ and $\psi _ { - } ( t , x ) = ( t - t _ { 0 } ) - | x - x _ { 0 } |$ , whose level surfaces $\psi _ { \pm } = 0$ correspond to forward and backward light cones with their vertex at $p = ( t _ { 0 } , x _ { 0 } )$ . These represent, physically, the union of all light rays emanating from a point source at p. The light rays are given by the equation $( t - t _ { 0 } ) \omega = ( x - x _ { 0 } )$ , for $\omega \in \mathbb { R } ^ { 3 }$ with $| \omega | = 1$ , and are precisely the (t, x) components of the bicharacteristic curves of the Hamilton–Jacobi equations (38). More generally, the characteristics of the linear wave equation

$$
a ^ {0 0} (t, x) \partial_ {t} ^ {2} u - \sum_ {i, j} a ^ {i j} (t, x) \partial_ {i} \partial_ {j} u = 0, \tag {39}
$$

with $a ^ { 0 0 } > 0$ and $\boldsymbol { a } ^ { i j }$ satisfying (35), are given by the Hamilton–Jacobi equations:

$$
- a ^ {0 0} (t, x) (\partial_ {t} \psi) ^ {2} + a ^ {i j} (x) \partial_ {i} \psi \partial_ {j} \psi = 0
$$

or, equivalently,

$$
\partial_ {t} \psi = \pm \left((a ^ {0 0}) ^ {- 1} \sum_ {i, j} a ^ {i j} (x) \partial_ {i} \psi \partial_ {j} \psi\right) ^ {1 / 2}. \tag {40}
$$

The bicharacteristics of the corresponding Hamiltonian systems are called bicharacteristic curves of (39).

Remark. In the case of the first-order scalar equations (17) we have seen how knowledge of characteristics can be used to find, implicitly, general solutions. We have also seen that singularities propagate only along characteristics. In the case of second-order equations the characteristics are not sufficient to solve the equations, but they continue to provide important information, such as how the singularities propagate. For example, in the case of the wave equation $\boldsymbol { \boxed { \ d } u } = 0$ with smooth initial data $u _ { 0 } , u _ { 1 }$ everywhere except at a point $p = ( t _ { 0 } , x _ { 0 } )$ , the solution u has singularities present at all points of the light cone $- ( t - t _ { 0 } ) ^ { 2 } + | x - x _ { 0 } | ^ { 2 } = 0$ with vertex at p. A more refined version of this fact shows that the singularities propagate along bicharacteristics. The general principle here is that singularities propagate along characteristic hypersurfaces of a PDE. Since this is a very important principle, it pays to give it a more precise formulation that extends to general boundary conditions, such as the Dirichlet condition for (1).

Propagation of singularities. If the boundary conditions or the coefficients of a PDE are singular at some point p, and otherwise smooth (or real analytic) everywhere in some small neighborhood V of p, then a solution of the equation cannot be singular in V except along a characteristic hypersurface passing through p. In particular, if there are no such characteristic hypersurfaces, then any solution of the equation must be smooth (or real analytic) at every point of V other than $p .$

Remarks. (i) The heuristic principle mentioned above is invalid, in general, at large scales. Indeed, as we have shown in the case of the Burgers equation, solutions to nonlinear evolution equations can develop new singularities whatever the smoothness of the initial conditions. Global versions of the principle can be formulated for linear equations based on the bicharacteristics of the equation. See (iii) below.

(ii) According to the principle, it follows that any solution of the equation $\Delta u \ = \ f ,$ , satisfying the boundary condition $u | _ { \partial D } = u _ { 0 }$ with a boundary value $u _ { 0 }$ that merely has to be continuous, is automatically smooth everywhere in the interior of D provided that f itself is smooth there. Moreover, the solution is real analytic if $f$ is real analytic.

(iii) More precise versions of this principle, which plays a fundamental role in the general theory, can be given for linear equations. In the case of the general wave equation (39), for example, one can show that singularities propagate along bicharacteristics. These are the bicharacteristic curves associated with the Hamilton– Jacobi equation (40).

# 2.4 The Cauchy–Kovalevskaya Theorem

In the case of ODEs we have seen that a noncharacteristic IVP always admits solutions locally (that is, in some time interval about a given point). Is there a higherdimensional analogue of this fact? The answer is yes, provided that we restrict ourselves to the real-analytic situation, which is covered by an appropriate extension of the Cauchy–Kovalevskaya theorem. More precisely, one can consider general quasilinear equations, or systems, with real-analytic coefficients, real-analytic hypersurfaces , and appropriate real-analytic initial data on .

Theorem (Cauchy–Kovalevskaya (CK)). If all the realanalyticity conditions made above are satisfied and if the initial hypersurface H is noncharacteristic at $x _ { 0 } , ^ { 4 }$ then in some neighborhood of $x _ { 0 }$ there is a unique real-analytic solution u(x) that satisfies the system of equations and the corresponding initial conditions.

In the special case of linear equations, an important companion theorem, due to Holmgren, asserts that the analytic solution given by the CK theorem is unique in the class of all smooth solutions and smooth noncharacteristic hypersurfaces H . The CK theorem shows that, given the noncharacteristic condition and the analyticity assumptions, the following straightforward way of finding solutions works: look for a formal expansion of the kind $\begin{array} { r } { \boldsymbol { u } ( \boldsymbol { x } ) = \sum _ { \alpha } C _ { \alpha } ( \boldsymbol { x } - \boldsymbol { x } _ { 0 } ) ^ { \alpha } } \end{array}$ by determining the constants $C _ { \alpha }$ recursively from simple algebraic formulas arising from the equation and initial conditions on . More precisely, the theorem ensures that the naive expansion obtained in this way converges in a small neighborhood of $x _ { 0 } \in \mathcal H$ .

It turns out, however, that the analyticity conditions required by the CK theorem are much too restrictive, and therefore the apparent generality of the result is misleading. A first limitation becomes immediately apparent when we consider the wave equation $\boldsymbol { \sqcup } \boldsymbol { u } = 0 ,$ . A fundamental feature of this equation is finite speed of propagation, which means, roughly speaking, that if at some time t a solution u is zero outside some bounded set, then the same must be true at all later times. However, analytic functions cannot have this property unless they are identically zero (see some fundamental mathematical definitions [I.3 §5.6]). Therefore, it is impossible to discuss the wave equation properly within the class of real-analytic solutions. A related problem, first pointed out by hadamard [VI.65], concerns the impossibility of solving the Cauchy problem, in many important cases, for arbitrary smooth nonanalytic data. Consider, for example, the Laplace equation $\Delta u = 0 \mathrm { i n } \mathbb { R } ^ { d }$ . As we have established above, any hypersurface  is noncharacteristic, yet the Cauchy problem $u | _ { \mathcal { H } } = u _ { 0 } , n [ u ] | _ { \mathcal { H } } = u _ { 1 }$ , for arbitrary smooth initial conditions $u _ { 0 } , u _ { 1 }$ , may admit no local solutions in a neighborhood of any point of ${ \mathcal { H } } .$ Indeed, take H to be the hyperplane $x _ { 1 } = 0$ and assume that the Cauchy problem can be solved for given nonanalytic smooth data in a domain that includes a closed ball B centered at the origin. The corresponding solution can also be interpreted as the solution to the Dirichlet problem in $B ,$ with the values of u prescribed on the boundary ∂B. But this, according to our heuristic principle (which can easily be made rigorous in this case), must be real analytic everywhere in the interior of $B ,$ contradicting our assumptions about the initial data.

On the other hand, the Cauchy problem for the wave equation $\boldsymbol { \Sigma } \boldsymbol { u } \ = \ \boldsymbol { 0 }$ in $\mathbb { R } ^ { d + 1 }$ has a unique solution for any smooth initial data $u _ { 0 } , \ u _ { 1 }$ that is prescribed on a spacelike hypersurface. This means a hypersurface $\psi ( t , x ) = 0$ such that at every point $p = ( t _ { 0 } , x _ { 0 } )$ that belongs to it the normal vector at p lies inside the light cone (either in the future direction or in the past direction). To say this analytically,

$$
\left| \partial_ {t} \psi (\boldsymbol {p}) \right| > \left(\sum_ {i = 1} ^ {d} \left| \partial_ {i} \psi (\boldsymbol {p}) \right| ^ {2}\right) ^ {1 / 2}. \tag {41}
$$

This condition is clearly satisfied by a hyperplane of the form $t \ = \ t _ { 0 }$ , but any other hypersurface close to this is also spacelike. By contrast, the IVP is ill-posed for a timelike hypersurface, i.e., a hypersurface for which

$$
\left| \partial_ {t} \psi (p) \right| <   \left(\sum_ {i = 1} ^ {d} \left| \partial_ {i} \psi (p) \right| ^ {2}\right) ^ {1 / 2}.
$$

That is, we cannot, for general non-real-analytic initial conditions, find a solution of the IVP. An example of a timelike hypersurface is given by the hyperplane $x ^ { 1 } =$ 0. Let us explain the term “ill-posed” more precisely.

Definition. A given problem for a PDE is said to be well-posed if both existence and uniqueness of solutions can be established for arbitrary data that belongs to a specified large space of functions, which includes the class of smooth functions.5 Moreover, the solutions must depend continuously on the data. A problem that is not well-posed is called ill-posed.

The continuous dependence on the data is very important. Indeed, the IVP would be of little use if very small changes in the initial conditions resulted in very large changes in the corresponding solutions.

# 2.5 Standard Classification

The different behavior of the Laplace and wave equations mentioned above illustrates the fundamental difference between ODEs and PDEs and the illusory generality of the CK theorem. Given that these two equations are so important in geometric and physical applications, it is of great interest to find the broadest classes of equations with which they share their main properties. The equations modeled by the Laplace equation are called elliptic, while those modeled by the wave equation are called hyperbolic. The other two important models are the heat equation (see (2)) and the Schrödinger equation (see (4)). The general classes of equations that they resemble are called parabolic and dispersive, respectively.

Elliptic equations are the most robust and the easiest to characterize: they are the ones that admit no characteristic hypersurfaces.

Definition. A linear, or quasilinear, $N \times N$ system with no characteristic hypersurfaces is called elliptic.

Equations of type (32) whose coefficients $a ^ { i j }$ satisfy condition (35) are clearly elliptic. The minimal surface equation (7) is also elliptic. It is also easy to verify that the Cauchy–Riemann system (12) is elliptic. As was pointed out by Hadamard, the IVP is not well-posed for elliptic equations. The natural way of parametrizing the set of solutions to an elliptic PDE is to prescribe conditions for u, and some of its derivatives (the number of derivatives will be roughly half the order of the equation) at the boundary of a domain $D \subset \mathbb { R } ^ { n }$ . These are called boundary-value problems (BVPs). A typical example is the Dirichlet boundary condition u|∂D = u0 for the Laplace equation $\Delta u \ = \ 0$ in a domain $D \subset \mathbb { R } ^ { n }$ . One can show that, if the domain D satisfies certain mild regularity assumptions and the boundary value $u _ { 0 }$ is continuous, then this problem admits a unique solution that depends continuously on $u _ { 0 } .$ . We say that the Dirichlet problem for the Laplace equation is wellposed. Another well-posed problem for the Laplace equation is given by the Neumann boundary condition $n [ u ] | _ { \partial D } = f ,$ , where n is the exterior unit normal to the boundary. This problem is well-posed for all continuous functions f defined on ∂D with zero mean average. A typical problem of general theory is to classify all well-posed BVPs for a given elliptic system.

As a consequence of our propagation-of-singularities principle, we deduce, heuristically at least, the following general fact:

Classical solutions of elliptic equations with smooth (or real-analytic) coefficients in a regular domain D are smooth (or real analytic) in the interior of D, whatever the degree of smoothness of the boundary conditions.6

Hyperbolic equations are, essentially, those for which the IVP is well-posed. In that sense, they provide the natural class of equations for which one can prove a result similar to the local existence theorem for ODEs. More precisely, for each sufficiently regular set of initial conditions there is a unique solution. We can thus think of the Cauchy problem as a natural way of parametrizing the set of all solutions to the equations.

The definition of hyperbolicity depends, however, on the particular hypersurface we are considering as the initial hypersurface. Thus, in the case of the wave equation $\boldsymbol { \sqcup } \boldsymbol { u } = 0 ,$ , the standard IVP

$$
u (0, x) = u _ {0} (x), \quad \partial_ {t} u (0, x) = u _ {1}
$$

is well-posed. This means that for any smooth initial data $u _ { 0 } , u _ { 1 }$ we can find a unique solution of the equation, which depends continuously on u0, u1. As we have already mentioned, the IVP for $\boxed { \ d } u = 0$ remains wellposed if we replace the initial hypersurface t = 0 by any spacelike hypersurface $\psi ( t , x ) = 0$ (see (41)). However, it fails to be well-posed for timelike hypersurfaces, for which there may be no solution with prescribed, nonanalytic, Cauchy data.

It is more difficult to give algebraic conditions for hyperbolicity. Roughly speaking, hyperbolic equations are at the opposite end of the spectrum from elliptic equations: whereas elliptic equations have no characteristic hypersurfaces, hyperbolic equations have as many as possible passing through any given point. One of the most useful classes of hyperbolic equations, which includes most of the important known examples, consists of equations of the form

$$
\begin{array}{r l} A ^ {0} (t, x, u) \partial_ {t} u + \sum_ {i = 1} ^ {d} A _ {i} (t, x, u) \partial_ {i} u & = F (t, x, u), \\ u | _ {\mathcal {H}} & = u _ {0}, \end{array} \tag {42}
$$

where all the coefficients $A ^ { 0 } , A ^ { 1 } , \ldots , A ^ { d }$ are symmetric $N \times N$ matrices and  is given by $\psi ( t , x ) = 0$ . Such a system is well-posed provided that the matrix

$$
A ^ {0} (t, x, u) \partial_ {t} \psi (t, x) + \sum_ {i = 1} ^ {d} A _ {i} (t, x, u) \partial_ {i} \psi (t, x) \tag {43}
$$

is positive definite. A system (42) that satisfies these conditions is called symmetric hyperbolic. In the particular case when $\psi ( t , x ) = t ,$ , the condition (43) becomes

$$
(A ^ {0} \xi , \xi) \geqslant c | \xi | ^ {2} \quad \forall \xi \in \mathbb {R} ^ {N}.
$$

The following is a fundamental result in the theory of general hyperbolic equations. It is called the local existence and uniqueness of solutions for symmetric hyperbolic systems.

Theorem (fundamental theorem for hyperbolic equations). The IVP (42) is locally well-posed for symmetric hyperbolic systems with sufficiently smooth $A , F ,$ and H and sufficiently smooth initial conditions $u _ { 0 } .$ . In other words, if the appropriate smoothness conditions are satisfied, then for any point $p \in { \mathcal { H } }$ there is a small neighborhood D of $p ^ { 7 }$ inside which there is a unique, continuously differentiable solution u.

Remarks. (i) The local character of the theorem is essential, just as it was for the general propagationof-singularities principle discussed earlier, since the result cannot be globalized in the particular case of the Burgers equation (21), which fits trivially into the framework of general nonlinear symmetric hyperbolic systems. A precise version of the theorem above gives a lower bound on how large D can be.

(ii) The proof of the theorem is based on a variation of the Picard iteration method that we encountered earlier for ODEs. One starts by taking $u _ { ( 0 ) } = u _ { 0 }$ in a neighborhood of H . Then one defines functions $u _ { ( n ) }$ recursively as follows:

$$
\begin{array}{l} A ^ {0} (t, x, u _ {(n - 1)}) \partial_ {t} u _ {(n)} + \sum_ {i = 1} ^ {d} A _ {i} (t, x, u _ {(n - 1)}) \partial_ {i} u _ {(n)} \\ = F (t, x, u _ {(n - 1)}), \quad u _ {(n)} | \mathcal {H} = u _ {0}. \\ \end{array}
$$

Notice that at each stage of the iteration we have to solve a linear equation. Linearization is an extremely important tool in studying nonlinear PDEs. We can almost never understand their behavior without linearizing them around important special solutions. Thus, almost invariably, hard problems in nonlinear PDEs reduce to understanding specific problems in linear PDEs.

(iii) To implement the Picard iteration method we need to get precise estimates concerning $u _ { ( n ) }$ in terms of $u _ { ( n - 1 ) }$ . This step requires energy type a priori estimates, which we will discuss in section 3.3.

Another important property of hyperbolic equations (which is not shared by elliptic, parabolic, or dispersive equations) is finite speed of propagation, which was mentioned earlier in the case of the wave equation (3). Consider this simple case again. The IVP can be solved explicitly by the so-called Kirchhoff formula. The formula allows us to conclude that if the initial data at $t \ = \ 0$ is zero outside a ball $B _ { a } ( x _ { 0 } )$ of radius $a > 0$ centered at $x _ { 0 } \in \mathbb { R } ^ { 3 }$ , then at time $t \ > \ 0$ the solution u is zero outside the ball $B _ { a + c t } ( x _ { 0 } )$ . In general, finite speed of propagation can best be formulated in terms of domains of dependence and influence of hyperbolic equations (see the online version for general definitions).

Hyperbolic PDEs play a fundamental role in physics, as they are intimately tied to the relativistic nature of the modern theory of fields. Equations (3), (5), (13) are the simplest examples of linear field theories, and they are manifestly hyperbolic. Other basic examples appear in gauge field theories such as maxwell’s equations [IV.13 §1.1] $\partial ^ { \alpha } F _ { \alpha \beta } = 0$ or the Yang–Mills equations $D ^ { \alpha } F _ { \alpha \beta } = 0$ . Finally, the Einstein equations (14) are also hyperbolic.8 Other important examples of hyperbolic equations arise in the physics of elasticity and inviscid fluids. As examples of the latter, the Burgers equation (21) and the compressible Euler equation are hyperbolic.

Elliptic equations, on the other hand, appear naturally in describing time-independent, or more generally steady-state, solutions of hyperbolic equations. Elliptic equations can also be derived, directly, by well-defined variational principles [III.94].

Finally, a few words about parabolic equations and Schrödinger-type equations, which are intermediate between the elliptic and hyperbolic ones. Large classes of useful equations of these types are given by

$$
\partial_ {t} u - L u = f \tag {44}
$$

and

$$
\mathrm{i} \partial_ {t} u + L u = f, \tag {45}
$$

respectively, where L is an elliptic second-order operator. One looks for solutions $u = u ( t , x )$ , defined for $t \geqslant t _ { 0 }$ , with the prescribed initial condition

$$
u (t _ {0}, x) = u _ {0} (x) \tag {46}
$$

on the hypersurface $t \ = \ t _ { 0 } .$ . Strictly speaking, this hypersurface is characteristic, since the order of the equation is 2 and we cannot determine $\partial _ { t } ^ { 2 }$ u at $t \ = \ t _ { 0 }$ directly from the equation. Yet this is not a serious problem; we can still determine $\partial _ { t } ^ { 2 }$ u formally by differentiating the equation with respect to $\partial _ { t } .$ . Thus, the IVP (44) (or (45)) with initial condition (46) is well-posed, but not quite in the same sense as for hyperbolic equations. For example, the heat equation $- \partial _ { t } u + \Delta u$ is well-posed for positive t but ill-posed for negative t. The heat equation may also not have unique solutions for the IVP unless we make assumptions about how fast the initial data is allowed to grow at infinity. One can also show that the characteristic hypersurfaces of the equation (44) are all of the form, and therefore parabolic equations are quite similar to elliptic equations. For example, one can show that if the coefficients $a ^ { i j }$ and $f$ are smooth (or real analytic), then the solution u must be smooth (or real analytic in x) for $t > t _ { 0 }$ even if the initial data $u _ { 0 }$ is not smooth, which is consistent with our propagation-of-singularities principle. The heat equation smooths out initial conditions. It is for this reason that the heat equation is useful in many applications. In physics, parabolic PDEs arise whenever diffusion or dissipation phenomena are important, while in geometry and calculus of variations, parabolic PDEs often arise as gradient flows of positive-definite functionals. Ricci flow (16) can also be viewed as a parabolic PDE, after a suitable change of coordinates.

Dispersive PDEs, of which the Schrödinger equation (4) is a fundamental example, are evolution equations that behave analogously to hyperbolic PDEs in many respects. For instance, the IVP tends to be locally well-posed both forward and backward in time. However, solutions to dispersive PDEs do not propagate along characteristic surfaces. Instead, they move at speeds that are determined by their spatial frequency; in general, high-frequency waves tend to propagate at much greater speeds than low-frequency waves, which eventually leads to a dispersion of the solution into increasingly large areas of space. In fact, the speed of propagation of solutions is typically infinite. This behavior also differs from that of parabolic equations, which tend to dissipate the high-frequency components of a solution (sending them to zero) rather than dispersing them. In physics, dispersive equations arise in quantum mechanics: they are the nonrelativistic limit c of relativistic equations and they are also approximations to model certain types of fluid behavior. For instance, the korteweg–de vries equation [III.49],

$$
\partial_ {t} u + \partial_ {x} ^ {3} u = 6 u \partial_ {x} u,
$$

is a dispersive PDE that models the behavior of smallamplitude waves in a shallow canal.

# 2.6 Special Topics for Linear Equations

The greatest successes of the general theory have been in connection with linear equations, especially those with constant coefficients, for which Fourier analysis provides an extremely powerful tool. While the related issues of classification, well-posedness, and propagation of singularities have dominated the study of linear equations, there are other issues of interest as well, including the following.

# 2.6.1 Local Solvability

This is the problem of determining the conditions on a linear operator  and given data f under which the equation (9) is locally solvable. The Cauchy–Kovalevskaya theorem gives a criterion for local solvability when f and the coefficients of  are real analytic, but it is a remarkable phenomenon that when one relaxes this assumption slightly, asking for $f$ to be smooth rather than real analytic, serious obstructions to local solvability appear. For instance, the Lewy operator

$$
\mathcal {P} [ u ] (t, z) = \frac {\partial u}{\partial \bar {z}} (t, z) - \mathrm{i} z \frac {\partial u}{\partial t} (t, z),
$$

defined on complex-valued functions $u : \mathbb { R } \times \mathbb { C } \to \mathbb { C } ,$ has the property that equation (9) is locally solvable for real-analytic f but not for “most” smooth $f .$ The Lewy operator is intimately connected to the tangential Cauchy–Riemann equations on the Heisenberg group in $\mathbb { C } ^ { 2 } .$ . It was discovered in the study of the restriction of the two-dimensional analogue of the Cauchy–Riemann operator $\mathcal { P }$ to a quadric in $\mathbb { C } ^ { 2 } .$ . This example was the starting point for the theory of local solvability, whose goal is to characterize linear equations that are locally solvable. The theory of Cauchy–Riemann manifolds— which has its origin in the study of restrictions of the Cauchy–Riemann equations (in higher dimensions) to real hypersurfaces, each of which comes with an associated “tangential Cauchy–Riemann complex”—is another extremely rich source of examples of interesting linear PDEs, which do not fit into the standard classification.

# 2.6.2 Unique Continuation

This concerns various ill-posed problems where solutions may not always exist, but one still has uniqueness. A fundamental example is that of analytic continuation: two holomorphic functions on a connected domain D that agree on a nondiscrete set (such as a disk or an interval) must necessarily agree everywhere on D. This fact can be viewed as a unique continuation result for the Cauchy–Riemann equations (12). Another example in a similar spirit is Holmgren’s theorem, which asserts that solutions to a linear PDE (9) that has realanalytic coefficients and data are unique, even in the class of smooth functions. More generally, the study of ill-posed problems (such as the wave equation with prescribed data on a timelike surface rather than a spacelike one) arises naturally in connection with control theory.

# 2.6.3 Spectral Theory

There is no way I can even begin to give an account of this theory, which is of fundamental importance not only to quantum mechanics and other physical theories, but also to geometry and analytic number theory [IV.2]. Just as a matrix A can often be analyzed through its eigenvalues and eigenvectors [I.3 §4.3] by the tools of linear algebra, one can learn much about a linear differential operator and its associated PDE by understanding that operator’s spectrum [III.86] and eigenfunctions with the help of tools from functional analysis [IV.15]. A typical problem in spectral theory is the eigenvalue problem in $\mathbb { R } ^ { d }$ :

$$
- \Delta u (x) + V (x) u (x) = \lambda u (x).
$$

A function u that is localized in space (for example, by being bounded in the $L ^ { 2 } ( \mathbb { R } ^ { d } )$ )-norm) and that satisfies this equation is mapped by the linear operator $- \Delta + V$ to the function λu: we say that u is an eigenfunction with eigenvalue λ.

Suppose that we have an eigenfunction u and let $\phi ( t , x ) ~ = ~ \mathrm { e } ^ { - \mathrm { i } \lambda t } u ( x )$ . It is easy to check that $\phi$ is a solution of the Schrödinger equation

$$
\mathrm{i} \partial_ {t} \phi + \Delta \phi - V \phi = 0. \tag {47}
$$

Moreover, it has a very special form. Such solutions are called bound states of the physical system described by (47). The eigenvalues λ, which form a discrete set, correspond to the quantum energy levels of the system. They are very sensitive to the choice of potential V . The inverse spectral problem is also important: can one determine the potential V from knowledge of the corresponding eigenvalues? The eigenvalue problem can be studied in considerable generality by replacing the operator $- \Delta + V$ with other elliptic operators. For instance, in geometry it is important to study the eigenvalue problem for the Laplace–Beltrami operator, which is the natural generalization of the Laplace operator from $\mathbb { R } ^ { n }$ to general riemannian manifolds [I.3 §6.10]. When the manifold has some arithmetic structure (for instance, if it is the quotient of the upper half-plane by a discrete arithmetic group), this problem is of major importance in number theory, leading, for instance, to the theory of Hecke–Maas forms. A famous problem in differential geometry (“can you hear the shape of a drum?”) is to characterize the metric on a compact surface from the spectral properties of the associated Laplace–Beltrami operator.

# 2.6.4 Scattering Theory

This theory formalizes the intuition from quantum mechanics that a potential which is small or localized is largely unable to “trap” a quantum particle, which is therefore likely to escape to infinity in a manner resembling that of a free particle. In the case of equation (47), solutions that scatter are those that behave freely as $t  \infty$ . That is, they behave like solutions to the free Schrödinger equation $\mathrm { i } \partial _ { t } \psi + \Delta \psi = 0 .$ A typical problem in scattering theory is to show that, if $V ( x )$ tends to zero sufficiently fast as $| x | \to \infty ,$ , all solutions, except the bound states, scatter as t .

# 2.7 Conclusions

In the analytic case, the CK theorem allows us to solve the IVP locally for very general classes of PDEs. We have a general theory of characteristic hypersurfaces of PDEs and a good general understanding of how they relate to propagation of singularities. We can also distinguish in considerable generality the fundamental classes of elliptic and hyperbolic equations and can define general parabolic and dispersive equations. The IVP for a large class of nonlinear hyperbolic systems can be solved locally in time, for sufficiently smooth initial conditions. Similar local-in-time results hold for general classes of nonlinear parabolic and dispersive equations. For linear equations a lot more can be done. We have satisfactory results concerning the regularity of solutions for elliptic and parabolic equations and a good understanding of the propagation of singularities for a large class of hyperbolic equations. Some aspects of spectral theory and scattering theory and problems of unique continuation can also be studied in considerable generality.

The main defect of the general theory concerns the passage from local to global. Important global features of special equations are too subtle to fit into a general scheme. Rather, each important PDE requires special treatment. This is particularly true for nonlinear equations: the long-term behavior of solutions is very sensitive to the special features of the equation at hand. Moreover, general points of view may obscure, through unnecessary technical complications, the main properties of the important special cases. A useful general framework is one that provides a simple and elegant treatment of a particular phenomenon, as is the case for symmetric hyperbolic systems and the phenomenon of local well-posedness and finite speed of propagation. However, it turns out that symmetric hyperbolic systems are simply too general for the study of more refined questions about the important examples of hyperbolic equations.

# 3 General Ideas

As one turns away from the general theory, one may be inclined to accept the pragmatic point of view described earlier, according to which PDEs is not a real subject but is rather a collection of subjects such as hydrodynamics, general relativity, several complex variables, elasticity, etc., each organized around a special equation. However, this rather widespread viewpoint has its own serious drawbacks. Even though specific equations have specific properties, the tools that are used to derive them are intimately related. In fact, there is an impressive body of knowledge relevant to all important equations, or at least large classes of them. Lack of space does not allow me to do anything more than enumerate them below.9

# 3.1 Well-Posedness

As is clear from the previous section, well-posed problems are at the heart of the modern theory of PDEs. Recall that these are problems that admit unique solutions for given smooth initial or boundary conditions, and that the corresponding solutions have to depend continuously on the data. It is this condition that leads to the classification of PDEs into elliptic, hyperbolic, parabolic, and dispersive equations. The first step in the study of a nonlinear evolution equation is a proof of a local-in-time existence and uniqueness theorem, similar to the one for ODEs. Ill-posedness, the counterpart of well-posedness, is also important in many applications. The Cauchy problem for the wave equation (3), with data on the timelike hypersurface z 0, is a typical example. Ill-posed problems appear naturally in control theory, as we have mentioned, and also in inverse scattering.

# 3.2 Explicit Representations and Fundamental Solutions

Our basic equations (2)–(5) can be solved explicitly. For example, the solution to the IVP for the heat equation in $\mathbb { R } _ { + } ^ { 1 + d }$ , that is, the problem of finding a function u that satisfies

$$
- \partial_ {t} u + \Delta u = 0, \quad u (0, x) = u _ {0} (x),
$$

for $t \geqslant 0 ,$ is given by

$$
u (t, x) = \int_ {\mathbb {R} ^ {d}} E _ {d} (t, x - y) u _ {0} (y) d y
$$

for a certain function $E _ { d } ,$ which is called the fundamental solution of the heat operator $- \partial _ { t } + \Delta$ . This function can be defined explicitly: when $t ~ \leqslant ~ 0$ it is $^ { 0 , }$ and when $t > 0$ it is given by the formula $E _ { d } ( t , x ) =$ (4πt) $^ { - d / 2 } \mathbf { e } ^ { - | x | ^ { 2 } / 4 t }$ . Observe that $E _ { d }$ satisfies the equation $( - \partial _ { t } + \Delta ) E = 0$ in both regions $t < 0$ and $t > 0 ,$ but it has a singularity at $t = 0 ,$ , which prevents it from satisfying the equation in the whole of $\mathbb { R } ^ { 1 + d }$ . In fact, we can check that for any function10 $\phi \in C _ { 0 } ^ { \infty } ( \mathbb { R } ^ { d + 1 } )$ , we have

$$
\int_ {\mathbb {R} ^ {d + 1}} E _ {d} (t, x) \left(\partial_ {t} \phi (t, x) + \Delta \phi (t, x)\right) \mathrm{d} t \mathrm{d} x = \phi (0, 0). \tag {48}
$$

In the language of distribution theory [III.18], formula (48) means that $E _ { d } ,$ as a distribution, satisfies the equation $( - \partial _ { t } + \Delta ) E _ { d } = \delta _ { 0 }$ , where $\delta _ { 0 }$ is the Dirac distribution in $\mathbb { R } ^ { 1 + d }$ supported at the origin. That is, $\delta _ { 0 } ( \phi ) = \phi ( 0 , 0 ) , \forall \phi \in C _ { 0 } ^ { \infty } ( \mathbb R ^ { d + 1 } )$ . A similar notion of fundamental solution can be defined for the Poisson, wave, Klein–Gordon, and Schrödinger equations.

A powerful method of solving linear PDEs with constant coefficients is based on the fourier transform [III.27]. For example, consider the heat equation $\partial _ { t } - \Delta u = 0$ in one space dimension, with initial condition $u ( 0 , x ) ~ = ~ u _ { 0 }$ . Define ${ \hat { u } } ( t , \xi )$ to be the Fourier transform of u relative to the space variable:

$$
\hat {u} (t, \xi) = \int_ {- \infty} ^ {+ \infty} \mathrm{e} ^ {- \mathrm{i} x \xi} u (t, x) \mathrm{d} x.
$$

It is easy to see that $\hat { u } ( t , \xi )$ satisfies the differential equation

$$
\partial_ {t} \hat {u} (t, \xi) = - \xi^ {2} \hat {u} (t, \xi), \quad \hat {u} (0, \xi) = \hat {u} _ {0} (\xi).
$$

This can be solved by a simple integration, which results in the formula $\hat { u } ( t , \xi ) = \hat { u } _ { 0 } ( \xi ) \mathrm { e } ^ { - t | \xi | ^ { 2 } }$ . Thus, with the help of the inverse Fourier transform, we derive a formula for u(t, x):

$$
u (t, x) = (2 \pi) ^ {- 1} \int_ {- \infty} ^ {+ \infty} \mathrm{e} ^ {\mathrm{i} x \xi} \mathrm{e} ^ {- t | \xi | ^ {2}} \hat {u} _ {0} (\xi) \mathrm{d} \xi .
$$

Similar formulas can be derived for our other basic evolution equations. For example, in the case of the wave equation $- \partial _ { t } ^ { 2 } u + \Delta u = 0$ in three dimensions, subject to the initial data $u ( 0 , x ) = u _ { 0 } , \hat { \sigma } _ { t } u ( 0 , x ) = 0 ,$ , we find that

$$
u (t, x) = (2 \pi) ^ {- 3} \int_ {\mathbb {R} ^ {3}} \mathrm{e} ^ {\mathrm{i} x \xi} \cos (t | \xi |) \hat {u} _ {0} (\xi) \mathrm{d} \xi . \tag {49}
$$

After some work, one can reexpress formula (49) in the form

$$
u (t, x) = \partial_ {t} \left((4 \pi t) ^ {- 1} \int_ {| x - y | = t} u _ {0} (y) \mathrm{d} a (y)\right), \tag {50}
$$

where da is the area element of the sphere x y t of radius t centered at x. This is the well-known Kirchhoff formula. By contrast with (49), the integration here is with respect to the physical variables t and x only. It is instructive to compare these two formulas. Using the Plancherel identity it is very easy to deduce from (49) the L2 bound

$$
\int_ {\mathbb {R} ^ {3}} | u (t, x) | ^ {2} \mathrm{d} x \leqslant C \| u _ {0} \| _ {L ^ {2} (\mathbb {R} ^ {3})} ^ {2},
$$

while the possibility of obtaining such a bound from (50) seems unlikely since the formula involves a derivative. On the other hand, (50) is perfect for giving us information about the domain of influence. Indeed, we can see immediately from the formula that if u0 is zero outside the ball $B _ { a } = \{ | x - x _ { 0 } | \leqslant a \}$ , then u(t, x) is zero outside the ball $B _ { a + | t | }$ for any time t. This fact does not seem at all transparent in the Fourier-based formula (49). The fact that different representations of solutions have different, even opposite, strengths and weaknesses has important consequences for constructing approximate solutions, or parametrices, for more complicated equations, such as linear equations with variable coefficients or nonlinear wave equations. There are two possible types of constructions: those in physical space, which mimic the physical-space formula (50), and those in Fourier space, which mimic the formula (49).

# 3.3 A Priori Estimates

Most equations cannot be solved explicitly. However, if we are interested in qualitative information about a solution, then it is not necessary to derive it from an exact formula. But how else, one might wonder, can we extract such information? A priori estimates are a very important technique for doing this.

The best-known examples are energy estimates, the maximum principle, and monotonicity arguments. The simplest example of the first type is the following identity (which is a very simple example of a so-called Bochner-type identity):

$$
\int_ {\mathbb {R} ^ {d}} | \partial^ {2} u (x) | ^ {2} d x = \int_ {\mathbb {R} ^ {d}} | \Delta u (x) | ^ {2} d x.
$$

The left-hand side is shorthand for

$$
\int_ {\mathbb {R} ^ {d}} \sum_ {1 \leqslant i, j \leqslant d} | \partial_ {i} \partial_ {j} u (x) | ^ {2} d x
$$

and the identity holds for all functions u that are twice continuously differentiable and tend to zero as |x| → ∞. This formula can be justified fairly simply by integrating by parts. As a consequence of the Bochner identity, we obtain the a priori estimate that if u is a smooth solution to the Poisson equation (6) with square-integrable data f , and if it tends to zero at infinity, then the square integral of its second derivatives is bounded:

$$
\int_ {\mathbb {R} ^ {d}} | \partial^ {2} u (x) | ^ {2} \mathrm{d} x \leqslant \int_ {\mathbb {R} ^ {d}} | f (x) | ^ {2} \mathrm{d} x <   \infty . \tag {51}
$$

Thus we obtain the qualitative fact that, on average (in a mean-square sense), u has “two more degrees of regularity” than f . 11 This is called an energy-type estimate because, in physical situations, the square of the L2-norm can often be interpreted as some type of kinetic energy.

The Bochner identity can be extended to more general Riemannian manifolds than Rd, although one then picks up some additional lower-order terms involving the curvature of those manifolds. Such identities play a major role in the theory of geometric PDEs on these manifolds.

Energy-type identities and estimates also exist for parabolic, dispersive, and hyperbolic PDEs. For instance, they play a fundamental role in demonstrating the local existence, uniqueness, and finite speed of propagation for hyperbolic PDEs with smooth initial data. Energy estimates become particularly powerful when combined with inequalities such as the Sobolev embedding inequality, which allows one to convert the $" L ^ { 2 } "$ information provided by these estimates into pointwise (or $^ { * } L ^ { \infty } )$ type information (see function spaces [III.29 §§2.4, 3]).

While energy identities and $L ^ { 2 }$ estimates (which, as in the above example, come from integration by parts) apply to all, or at least major classes of, PDEs, the maximum principle can be applied only to elliptic and parabolic PDEs. The following theorem is the simplest manifestation of it. Note that the theorem provides us with important quantitative information about solutions to the Laplace equation even in the absence of any explicit representation for them.

Theorem (maximum principle). Assume that u is a solution to the Laplace equation (1) on a bounded connected domain $D \in \mathbb { R } ^ { d }$ with a smooth boundary ∂D. Assume also that u is continuous on the closure of D and has continuous first and second partial derivatives in the interior of D. Then u must achieve its maximum and minimum values on the boundary. Moreover, if the maximum or minimum is also achieved at an interior point of D, then u must be constant in D.

The method is very robust and can easily be extended to a large class of second-order elliptic equations. It can also be extended to parabolic equations and systems, and plays a crucial role ${ \mathrm { i n } } ,$ for example, the study of Ricci flow.

Let us briefly mention some other important classes of a priori estimates. The Sobolev inequalities, which are of prime importance in elliptic equations, have several counterparts in linear and nonlinear hyperbolic and dispersive equations, including the Strichartz estimates and bilinear estimates. In connection with ill-posed problems and unique continuation, Carleman estimates play a fundamental role. Finally, several a priori estimates arising from monotonicity formulas12—such as virial identities, Pohozaev identities, or Morawetz inequalities—can be used to establish the breakdown of regularity or the blow-up of solutions to some nonlinear equations, and to guarantee global existence and decay of solutions to others.

To summarize, it is not much of an exaggeration to say that a priori estimates play a fundamental role in more or less every aspect of the modern theory of PDEs.

# 3.4 Bootstrap and Continuity Arguments

The bootstrap argument is a method, or rather a powerful general philosophy, to derive a priori estimates for nonlinear equations. According to this philosophy we start by making educated assumptions about the solutions we are trying to describe. These assumptions allow us to think of the original nonlinear problem as a linear one whose coefficients satisfy properties consistent with the assumptions. We may then use linear methods, based on other a priori estimates that we already know, to try to show that the solutions to this linear problem behave as well as we have postulated— in fact, even better. One can characterize this powerful method, which allows us to use linear theory without actually having to linearize the equation, as a conceptual linearization. It can also be regarded as a continuity argument relative to some parameter, which might be the natural time parameter of an evolution problem, but it could also be an artificial parameter which we have the freedom to introduce ourselves. This latter situation is typical of applications to nonlinear elliptic equations. In the online version of this article we provide a few examples to illustrate the method in both cases.

# 3.5 The Method of Generalized Solutions

Since a PDE involves differentiation, it might seem obvious that in any discussion of PDEs we should restrict our attention to differentiable functions. However, it is possible to generalize the notion of differentiation so that it makes sense for a wider class of functions, and even for function-like objects, such as distributions, that are not functions at all. This allows us to make sense of a PDE in a broader context, and admits the possibility of generalized solutions.

The best way to introduce generalized solutions in PDEs and explain why they are important is through the Dirichlet principle. This originates in the observation that, out of all functions that are defined on a bounded domain $D ~ \subset ~ \mathbb { R } ^ { d }$ , that satisfy prescribed Dirichlet boundary condition u $\partial D = f$ , and that live in an appropriate functional space X, the functions u that minimize the Dirichlet integral (or Dirichlet functional)

$$
\| u \| _ {D r} ^ {2} = \frac {1}{2} \int_ {D} | \nabla u | ^ {2} = \frac {1}{2} \sum_ {i = 1} ^ {d} \int_ {D} | \partial_ {i} u | ^ {2} \tag {52}
$$

are the harmonic functions (that is, solutions of the equation $\Delta u = 0 )$ . It was riemann [VI.49] who first had the idea of trying to use this fact to solve the Dirichlet problem: in order to find a solution u to the problem

$$
\Delta u = 0, \quad u | _ {\partial D} = u _ {0}, \tag {53}
$$

one should find (by some means other than solving the Dirichlet problem) a function u that minimizes the Dirichlet integral while equaling $u _ { 0 }$ on $\partial D .$ To do this, one must specify the set by functions, or rather the function space, over which the minimization is taking place. The history of how this choice was made is a fascinating one. A natural choice is $X = C ^ { 1 } ( \bar { D } )$ , the space of continuously differentiable functions on D¯, where the norm of a function v is

$$
\left\| \boldsymbol {v} \right\| _ {C ^ {1} (\bar {D})} = \sup _ {x \in D} (| \boldsymbol {v} (x) | + | \partial \boldsymbol {v} (x) |).
$$

In particular, the Dirichlet norm $\| \boldsymbol { \nu } \| _ { D r }$ is finite when v belongs to this space. In fact, Riemann chose $X \ =$ $C ^ { 2 } ( { \bar { D } } )$ (a similar space but designed for twice continuously differentiable functions). This bold but flawed attempt was followed by a penetrating criticism by weierstrass [VI.44], who showed that the functional does not have to achieve its minimum in either $C ^ { 2 } ( \bar { D } )$ or $C ^ { 1 } ( \bar { D } )$ . However, Riemann’s basic idea was revived, and it eventually triumphed after a long and inspiring process that involved defining appropriate function spaces, introducing the notion of generalized solutions, and developing a regularity theory for them. (The precise formulation of the Dirichlet principle also requires the definition of sobolev spaces [III.29 §2.4].)

Let us briefly summarize the method, which has since been vastly extended so that it can be applied to a large class of $\mathrm { l i n e a r ^ { 1 3 } }$ and nonlinear elliptic and parabolic equations. It is based on two steps. In the first step one applies a minimization procedure. Although, as Weierstrass discovered, the natural function spaces may not contain functions that achieve the minimum, one can use such a procedure to find a generalized solution instead. This may not seem very interesting, since we were looking for a function that solves the Dirichlet problem (or one of the other problems to which the method can be applied). But this is where the second step comes in: it is sometimes possible to show that the generalized solution must in fact be a classical solution (that is, an appropriately smooth function) after all. This is the “regularity theory” mentioned earlier. In some situations, however, the generalized solution may turn out to have singularities and therefore not be regular. Then the challenge is to understand the nature of these singularities and to prove realistic partial regularity results. For instance, it is sometimes possible to prove that the generalized solution is smooth everywhere apart from in a small “exceptional set.”

Though generalized solutions are at their most effective for elliptic problems, their range of applicability encompasses all PDEs. For example, we have already seen that the fundamental solutions to the basic linear equations have to be interpreted as distributions, which are examples of generalized solutions.

The notion of generalized solutions has also proved successful for nonlinear evolution problems, such as systems of conservation laws in one space dimension. An excellent example is provided by the Burgers equation (21). As we have seen, solutions to $\partial _ { t } u + u \partial _ { x } u =$ 0 develop singularities in finite time no matter how smooth the initial conditions are. It is natural to ask whether solutions continue to make sense, as generalized solutions, even beyond the time when these singularities form. A natural notion of generalized solution is a function u such that

$$
\int_ {\mathbb {R} ^ {1 + 1}} \left(\partial_ {t} u + u \partial_ {x} u\right) \phi = 0
$$

for every smooth function $\phi$ that is zero outside a bounded set, since one can make sense of the integral even when u is not a differentiable function. Integrating this by parts (the first term with respect to t and the second with respect to x) one obtains the following formulation:

$$
\int_ {\mathbb {R} ^ {1 + 1}} u \partial_ {t} \phi + \frac {1}{2} \int_ {\mathbb {R} ^ {1 + 1}} u ^ {2} \partial_ {x} \phi = 0 \quad \forall \phi \in C _ {0} ^ {\infty} (\mathbb {R} ^ {1 + 1}).
$$

It can be shown that, under additional conditions called entropy conditions, the IVP for the Burgers equation admits a unique generalized solution that is global: that is, valid for every $t \in \mathbb { R }$ . Today we have a satisfactory theory of global solutions to a large class of hyperbolic systems of one-dimensional “conservation laws.” These systems, for which the above-mentioned theory applies, are called strictly hyperbolic.

For more complicated nonlinear evolution equations, the question of what constitutes a good concept of a generalized solution, though fundamental, is far murkier. For higher-dimensional evolution equations the first concept of a weak solution was introduced by Leray. Let us call a generalized solution weak if one cannot prove any type of uniqueness for it. This unsatisfactory situation may be temporary, i.e., the result of our technical inabilities, or unavoidable, in the sense that the concept itself is flawed. Leray was able to produce, by a compactness method, a weak solution of the IVP for the navier–stokes equations [III.23]. The great advantage of the compactness method (and its modern extensions, which can, in some cases, cleverly circumvent lack of compactness) is that it produces global solutions for all data. This is particularly important for supercritical or critical nonlinear evolution equations, which we will discuss later. For these we expect classical solutions to develop singularities in a finite time. The problem, however, is that one has very little control over such solutions. In particular, we do not know how to prove their uniqueness.14 Similar types of solutions were later introduced for other important nonlinear evolution equations. In most of the interesting cases of supercritical evolution equations, such as the Navier–Stokes equations, the usefulness of the types of weak solutions discovered so far remains undecided.

# 3.6 Microlocal Analysis, Parametrices, and Paradifferential Calculus

One of the fundamental difficulties of hyperbolic and dispersive equations is the interplay between geometric properties, which concern the physical space, and other properties, intimately tied to oscillations, that are best seen in Fourier space. Microlocal analysis is a general still-developing philosophy according to which one isolates the main difficulties by careful localizations in physical space or Fourier space or both. An important application of this point of view is the construction of parametrices for linear hyperbolic equations and their use in proving results about the propagation of singularities. Parametrices, as we have already mentioned, are approximate solutions of linear equations with variable coefficients, with error terms that are smoother. The paradifferential calculus is an extension of microlocal analysis to nonlinear equations. It allows one to manipulate the form of a nonlinear equation by taking account of how large and small frequencies interact, and it has achieved a remarkable technical versatility.

# 3.7 Scaling Properties of Nonlinear Equations

A PDE is said to have a scaling property if, whenever one rescales a solution in an appropriate way, one obtains another solution. Essentially, all basic nonlinear equations have well-defined scaling properties. Take, for example, the Burgers equation (21), $\partial _ { t } u + u \partial _ { x } u = 0 .$ If u is a solution of this equation, then so is the function $u _ { \lambda }$ defined by $u _ { \lambda } ( t , x ) = u ( \lambda t , \lambda x )$ . Similarly, if u is a solution of the cubic nonlinear Schrödinger equation in $\mathbb { R } ^ { d }$ ,

$$
\mathrm{i} \partial_ {t} u + \Delta u + c | u | ^ {2} u = 0, \tag {54}
$$

then so is $u _ { \lambda } ( t , x ) ~ = ~ \lambda u ( \lambda ^ { 2 } t , \lambda x )$ . The relationship between the nonlinear scaling of the equation and the a priori estimates available for solutions to the equations leads to an extremely useful classification of equations into subcritical, critical, and supercritical equations. This will be discussed in more detail in the next section. For the moment it suffices to say that subcritical equations are those for which the nonlinearity can be controlled by the existing a priori estimates of the equation, while supercritical equations are those for which the nonlinearity appears to be stronger. Critical equations are borderline. The definition of criticality and its relationship with the issue of regularity play a very important heuristic role in nonlinear PDEs. One expects supercritical equations to develop singularities and subcritical equations not to.

# 4 The Main Equations

In the previous section we argued that, while there is no hope of finding a general theory of all PDEs, there is nevertheless a wealth of general ideas and techniques that are relevant to the study of almost all important equations. In this section we indicate how it may be possible to identify the features that characterize the equations we call important.

Most of our basic PDEs can be derived from simple geometric principles, which happen to coincide with some of the underlying geometric principles of modern physics. These simple principles provide a unifying framework15 for the subject and help endow it with a sense of purpose and cohesion. They also explain why a very small number of linear differential operators, such as the Laplacian and the d’Alembertian, are all-pervasive.

Let us begin with the operators. The Laplacian is the simplest differential operator that is invariant under rigid motions of Euclidean space—a fact that we noted at the beginning of this article. This is important mathematically and physically: mathematically because it results in many symmetry properties and physically because many physical laws are themselves invariant under rigid motions. The d’Alembertian is, similarly, the simplest differential operator that is invariant under the natural symmetries, or Poincaré transformations, of Minkowski space.

Now let us turn to the equations. From the point of view of physics, the heat equation is basic because it is the simplest paradigm for diffusive phenomena, while the Schrödinger equation can be viewed as the Newtonian limit of the Klein–Gordon equation. The geometric framework of the former is Galilean space, which itself is simply the Newtonian limit of Minkowski space.16

From a mathematical point of view, the heat, Schrö- dinger, and wave equations are basic because the corresponding differential operators $\partial _ { t } - \Delta , ( 1 / \mathrm { i } ) \partial _ { t } - \Delta ,$ , and $\partial _ { t } ^ { 2 } - \Delta$ are the simplest evolution operators that can be built out of Δ. The wave operator, as just discussed, is basic in a deeper way because of the association between $\bigtriangledown = - \partial _ { t } ^ { 2 } + \Delta$ and the geometry of Minkowski space $\mathbb { R } ^ { 1 + n }$ . As for Laplace’s equation, one can view solutions to $\Delta \phi = 0$ as special time-independent solutions to $\bigstar \phi = 0$ . Appropriate invariant and local definitions of square roots of Δ and $\begin{array} { r } { \bigsqcup , \ : \mathrm { o r } \ : \bigsqcup - k ^ { 2 } } \end{array}$ , corresponding to “spinorial representations” of the Lorentz group, lead to the associated Dirac operators (see (13)). In the same vein we can associate with every Riemannian or Lorentzian manifold the operator $\Delta _ { g }$ or $\sqsubseteq _ { g } ,$ , respectively, or the corresponding Dirac operators. These equations inherit in a straightforward way the symmetries of the spaces on which they are defined.

# 4.1 Variational Equations

There is a general and extremely effective method for generating equations with prescribed symmetries that plays a fundamental role in both physics and geometry. One starts with a scalar quantity, called a Lagrangian, such as

$$
\mathcal {L} [ \phi ] = \sum_ {\mu , \nu = 0} ^ {3} m ^ {\mu \nu} \partial_ {\mu} \phi \partial_ {\nu} \phi - V (\phi), \tag {55}
$$

with $\phi$ a real-valued function defined on $\mathbb { R } ^ { 1 + 3 }$ and V some real function of $\phi$ such as, for example, $V ( \phi ) =$ $\phi ^ { 3 }$ . Here $\partial _ { \mu }$ denotes the partial derivatives with respect to the coordinates $x ^ { \mu } , \mu = 0 , 1 , 2 , 3$ , and $m ^ { \mu \nu } = m _ { \mu \nu }$ , as earlier, denotes the 4 × 4 diagonal matrix with diagonal entries $( - 1 , 1 , 1 , 1 )$ , associated with the Minkowski metric. We associate with ${ \mathcal { L } } [ \phi ]$ the so-called action integral:

$$
S [ \phi ] = \int_ {\mathbb {R} ^ {3 + 1}} \mathcal {L} [ \phi ].
$$

Notice that both ${ \mathcal { L } } [ \phi ]$ and $S [ \phi ]$ are invariant under translations and Lorentz transformations. In other words, if $T : \mathbb { R } ^ { 1 + 3 } \  \ \mathbb { R } ^ { 1 + 3 }$ is a function that does not change the metric and we define a new function by $\psi ( t , x ) = \phi ( T ( t , x ) )$ , then ${ \mathcal { L } } [ \phi ] = { \mathcal { L } } [ \psi ]$ and $S [ \phi ] =$ $S [ \psi ]$ .

We shall consider a function φ that minimizes the action integral. From this we wish to deduce that the derivative of  at $\phi ,$ , in some appropriate sense, is zero, and hence to deduce other properties about $\phi .$ But $\phi$ is a function that lives in an infinite-dimensional space, so we cannot talk about derivatives in a completely straightforward way. To deal with this problem, we define a compact variation of φ to be a smooth oneparameter family of functions $\phi ^ { ( s ) } : \mathbb { R } ^ { 1 + 3 }  \mathbb { R } ,$ , defined for each s in some interval $( - \epsilon , \epsilon )$ , such that $\phi ^ { ( 0 ) } ( x ) =$ $\phi ( x )$ for every $\textbf { \textit { x } } \in { \mathrm { ~ \mathbb ~ R ^ { 3 } ~ } }$ and $\phi ^ { ( s ) } ( x ) = \phi ( x )$ for every $( s , x )$ outside some bounded subset of $\mathbb { R } ^ { 1 + 3 }$ . This allows us to differentiate with respect to s.

Given such a variation, we denote the derivative $\mathrm { d } \phi ^ { ( s ) } / \mathrm { d } s | _ { s = 0 }$ by $\dot { \phi } .$

Definition. A field $\phi$ is said to be stationary with respect to ${ \mathcal { S } } { \mathrm { ~ i f } } ,$ for any compact variation $\phi ^ { ( s ) }$ of $\phi ,$ we have

$$
\left. \frac {\mathrm{d}}{\mathrm{d} s} S [ \phi^ {(s)} ] \right| _ {s = 0} = 0.
$$

The variational principle. The variational principle, or principle of least action, states that an acceptable solution of a given physical system must be stationary with respect to the action integral associated with the Lagrangian of the system.

The variational principle enables us to associate with the given Lagrangian a system of PDEs, obtained from the fact that φ is stationary, called the Euler–Lagrange equations. We illustrate this by showing that the nonlinear wave equation in $\mathbb { R } ^ { 1 + 3 }$ , namely

$$
\Box \phi - V ^ {\prime} (\phi) = 0, \tag {56}
$$

is the Euler–Lagrange equation associated with the Lagrangian (55). Given a compact variation $\phi ^ { ( s ) }$ of $\phi ,$ we set $S ( s ) = S [ \phi ^ { ( s ) } ]$ . Integration by parts gives

$$
\begin{array}{l} \left. \frac {\mathrm{d}}{\mathrm{d} s} S (s) \right| _ {s = 0} = \int_ {\mathbb {R} ^ {3 + 1}} \left[ - m ^ {\mu \nu} \partial_ {\mu} \dot {\phi} \partial_ {\nu} \phi - V ^ {\prime} (\phi) \dot {\phi} \right] \\ = \int_ {\mathbb {R} ^ {3 + 1}} \dot {\phi} [ \square \phi - V ^ {\prime} (\phi) ]. \\ \end{array}
$$

In view of the action principle and the arbitrariness of φ˙ we infer that φ must satisfy equation (56). Thus (56) is indeed the Euler–Lagrange equation associated with the Lagrangian ${ \mathcal { L } } [ \phi ] = m ^ { \mu \nu } \partial _ { \mu } \phi \partial _ { \nu } \phi - V ( \phi )$ .

One can similarly show that the Maxwell equations of electromagnetism—along with their beautiful extensions to the Yang–Mills equations, wave maps, and the Einstein equations of general relativity—are also variational. That is, they too can be derived from a Lagrangian.

Remark. The variational principle asserts only that the acceptable solutions of a given system are stationary: in general, we have no reason to expect that the desired solutions minimize or maximize the action integral. Indeed, this fails to be the case for systems that have a time dependence, such as the Maxwell equations, Yang–Mills equations, wave maps, and Einstein equations.

However, there is a large class of variational problems, corresponding to time-independent physical systems or geometric problems, for which the desired solutions do turn out to be extremal. The simplest example is that of geodesics in a Riemannian manifold M, which are minimizers17 with respect to length. More precisely, the length functional takes a curve γ that passes through two fixed points of M and associates with it its length $L ( \gamma )$ , which plays the role of an action integral. In this case a geodesic is not just a stationary point for the functional but a minimum. We also saw earlier that, according to the Dirichlet principle, solutions to the Dirichlet problem (53) minimize the Dirichlet integral (52). Another example is provided by the minimal surface equation (7), the solutions of which are minimizers of the area integral.

The study of minimizers of various functionals, i.e., action integrals, is a venerable subject in mathematics that goes under the name of calculus of variations (see variational methods [III.94] for further discussion).

Associated with the variational principle is another fundamental principle. A conservation law for an evolution PDE is a law that says that some quantity, typically an integral quantity depending on the solution, must remain constant over time, for every solution of the equation.

Noether’s principle. To any continuous one-parameter group of symmetries of the Lagrangian there corresponds a conservation law for the associated Euler– Lagrange PDE.

Examples of such conservation laws are the familiar laws of conservation of energy, conservation of momentum, and conservation of angular momentum, all of which have important physical meaning. (The oneparameter group of symmetries for energy, for example, is just translations in time.) In the case of equation (56), the law of conservation of energy takes the form

$$
E (t) = E (0), \tag {57}
$$

where the quantity E(t), which equals

$$
\int_ {\Sigma_ {t}} \left(\frac {1}{2} (\partial_ {t} \phi) ^ {2} + \frac {1}{2} \sum_ {i = 1} ^ {3} (\partial_ {i} \phi) ^ {2} + V (\phi)\right) d x, \tag {58}
$$

is called the total energy at time t. (We write $\Sigma _ { t }$ for the set of all points $( t , x , y , z )$ as (x, y, z) ranges over $\mathbb { R } ^ { 3 } .$ ) Observe that (57) provides an extremely important a priori estimate for solutions to (56) in the case when $V \geqslant 0$ . Indeed, if the energy of the initial data at $t = 0$ is finite (that is, if $E ( 0 ) < \infty )$ , then

$$
\int_ {\Sigma_ {t}} \left(\frac {1}{2} (\partial_ {t} \phi) ^ {2} + \frac {1}{2} \sum_ {i = 1} ^ {3} (\partial_ {i} \phi) ^ {2}\right) d x \leqslant E (0).
$$

We say that the energy identity (57) is coercive, which means that it leads to an absolute bound on all solutions with finite initial energy.

# 4.2 The Issue of Criticality

For the most basic evolution equations of mathematical physics, there are typically no better a priori estimates known than those provided by the energy. Taking into account the scaling properties of the corresponding equations as well, one is led to the very important classification of our basic equations, mentioned earlier, into subcritical, critical, and supercritical equations. To see how this is done, consider again the nonlinear scalar equation φ − $V ^ { \prime } ( \phi ) \ = \ 0$ , and take $V ( \phi )$ to be $( 1 / ( p + 1 ) ) | \phi | ^ { p + 1 }$ . Recall that the energy integral is given by (58). If we assign to the spacetime variables the dimension of length, L, then the spacetime derivatives have dimension $L ^ { - 1 }$ and therefore  has the dimension of $L ^ { - 2 }$ . To be able to balance the left- and right-hand sides of the equation $\square \phi = | \phi | ^ { p - 1 } \phi$ , we need to assign a length scale to φ; we find this to be $L ^ { 2 / ( 1 - p ) }$ . Thus the energy integral,

$$
E (t) = \int_ {\mathbb {R} ^ {d}} \left(2 ^ {- 1} | \partial \phi | ^ {2} + | \phi | ^ {p + 1}\right) d x,
$$

has the dimension $L ^ { c } , c = d \mathrm { - } 2 + ( 4 / ( 1 - p ) )$ , with d corresponding to the volume element $\mathrm { d } x = \mathrm { d } x ^ { \mathrm { 1 } } \mathrm { d } x ^ { \mathrm { 2 } }$ $\mathrm { d } x ^ { d }$ , which scales like $L ^ { d } .$ . We say that the equation is subcritical if $c < 0 .$ , critical if $c = 0 ,$ , and supercritical if $c > 0 ,$ . Thus, for example, $\begin{array} { r } { \square \phi - \phi ^ { 5 } = 0 } \end{array}$ is critical in dimension $d = 3 ,$ . The same sort of dimensional analysis can be done for all our other basic equations. An evolutionary PDE is said to be regular if all smooth finite-energy initial conditions lead to global smooth solutions. It is conjectured that all subcritical equations are regular, but one expects supercritical equations to develop singularities. Critical equations are important borderline cases. The heuristic reason for this is that the nonlinearity tends to produce singularities while the coercive estimates prevent it. In subcritical equations the coercive estimates are stronger, while for supercritical equations it is the nonlinearity that is stronger. However, there may be other, more subtle a priori estimates that are not accounted for by our crude heuristic argument. Thus, some supercritical equations, such as the Navier–Stokes equations, may still be regular.

# 4.3 Other Equations

Many other familiar equations can be derived from the variational ones described above by the following procedures.

# 4.3.1 Symmetry Reductions

Sometimes a PDE is very hard to solve but becomes much easier if one places additional symmetry constraints on solutions. For example, if the PDE is rotation invariant and we look just for rotation-invariant solutions $u ( t , x )$ , then we can regard these solutions as functions of t and $r = \left| \boldsymbol { x } \right|$ , effectively reducing the dimension of the problem. By this procedure of symmetry reduction one can then derive a new PDE that is much simpler than the original one. Another, somewhat more general, way of obtaining simpler equations is to look for solutions that satisfy some further property. For instance, one can assume that they are stationary (that is, that they do not depend on the time variable), spherically symmetric, self-similar (which means that u(t, x) depends only on $x / t ^ { a } )$ , or traveling waves (which means that $u ( t , x )$ depends only on $x - \upsilon t$ for some fixed velocity vector v). Typically, the equations obtained by such reductions have a variational structure themselves. In fact, the symmetry reduction can be applied directly to the original Lagrangian.

# 4.3.2 The Newtonian Approximation and Other Limits

We can derive a large class of new equations as limits of the basic ones described above by taking one or more characteristic speeds to infinity. The most important example is the Newtonian limit, which is obtained by letting the velocity of light go to infinity. As we have already mentioned, the Schrödinger equation can be derived in this way from the linear Klein–Gordon equation. Similarly, we can derive the Lagrangians for the equations of nonrelativistic elasticity, fluid dynamics, or magnetohydrodynamics. It is an interesting fact that the nonrelativistic equations tend to look more messy than the relativistic ones. The simple geometric structure of the original equations gets lost in the limit. The remarkable simplicity of the relativistic equations is a powerful example of the importance of relativity as a unifying principle.

Once we are in the familiar world of Newtonian physics we can perform other well-known limiting procedures. The famous incompressible euler equations [III.23] are obtained by taking the limit of the general nonrelativistic fluid equations as the speed of sound tends to infinity. Various other limits are obtained relative to other characteristic speeds of the system or in connection with specific boundary conditions, such as the boundary-layer approximation in fluids. For example, in the limit as all characteristic speeds tend to infinity, the equations of elasticity turn into the familiar equations of a rigid body in classical mechanics.

# 4.3.3 Phenomenological Assumptions

Even after taking various limits and making symmetry reductions, the equations may still remain intractable. However, in various applications it makes sense to assume that certain quantities are sufficiently small to be neglected. This leads to simplified equations that could be called phenomenological18 in the sense that they are not derived from first principles.

Phenomenological equations are “toy equations” that are used to illustrate and isolate important physical phenomena in complicated systems. A typical way of generating interesting phenomenological equations is to try to write down the simplest model equation that still exhibits a particular feature of the original system. For instance, the self-focusing plane-wave effects of compressible fluids or elasticity can be illustrated by the simple-minded Burgers equation $u _ { t } + u u _ { x } ~ =$ 0. Nonlinear dispersive phenomena, typical of fluids, can be illustrated by the famous Korteweg–de Vries equation $u _ { t } + u u _ { x } + u _ { x x x } = 0$ . The nonlinear Schrö- dinger equation (54) provides a good model problem for nonlinear dispersive effects in optics.

If it is well chosen, a model equation can lead to basic insights into the original equation itself. For this reason, simplified model problems are also essential in the day-to-day work of the rigorous researcher into PDEs, who tests ideas on carefully selected model problems. It is crucial to emphasize that good results concerning the basic physical equations are rare; a very large percentage of important rigorous work in PDEs deals with simplified equations selected, for technical reasons, to isolate and focus our attention on some specific difficulties present in the basic equations.

In the above discussion we have not mentioned diffusive equations19 such as the Navier–Stokes equations. These are in fact not variational, and therefore do not quite fit into the above description. Though they could be viewed as phenomenological equations, they can also be derived from basic microscopic laws such as those governing the Newtonian–mechanical interactions of a very large number of particles N. In principle,20 the equations of continuum mechanics, such as the Navier–Stokes equations, could be derived by letting the number of particles N → ∞.

Diffusive equations also turn out to be very useful in connection with geometric problems. Geometric flows such as mean curvature, inverse mean curvature, harmonic maps, Gauss curvature, and Ricci flow are some of the best-known examples. Diffusive equations can often be interpreted as the gradient flow for an associated elliptic variational problem. They can be used to construct nontrivial stationary solutions to the corresponding stationary systems, in the limit as t  , or to produce foliations with remarkable properties, such as one that was used recently in the proof of a famous conjecture of Penrose. As we have already mentioned, this idea has recently found an extraordinary application in the work of Perelman, who has used Ricci flow to settle the three-dimensional Poincaré conjecture. One of his main new ideas was to interpret Ricci flow as a gradient flow.

# 4.4 Regularity or Breakdown

An additional source of unity for the subject of PDEs is the central role played by the problem of regularity or breakdown of solutions to the basic equations. It is intimately tied to the fundamental mathematical question of understanding what we actually mean by solutions and, from a physical point of view, to the issue of understanding the limits of validity of the corresponding physical theories. Thus, in the case of the Burgers equation, for example, the problem of singularities can be tackled by extending our concept of solutions to accommodate shock waves, which are solutions that are discontinuous across certain curves in the (t, x)-space. In this case one can define a function space of generalized solutions in which the IVP has unique, global solutions. Though the situation for more realistic physical systems is far less clear and far from being satisfactorily solved, the generally held opinion is that shockwave-type singularities can be accommodated without breaking the boundaries of the physical theory at hand. The situation for singularities in general relativity is radically different. The singularities one expects there are such that no continuation of solutions is possible without altering the physical theory itself. The prevailing opinion here is that only a gravitational quantum field theory could achieve this.

# 5 General Conclusions

What, then, is the modern theory of PDEs? As a first approximation, one could say that it is the pursuit of the following main goals.

(i) Understand the problem of evolution for the basic equations of mathematical physics. The most pressing issue in this regard is to understand when and how the local21 (with respect to time) smooth solutions of the basic equations develop singularities. A simple-minded criterion for distinguishing between regular theories and those that may admit singular solutions is given by the distinction between subcritical and supercritical equations. As mentioned earlier, it is widely believed that subcritical equations are regular and that supercritical equations are not. Indeed, many subcritical equations have been proved to be regular even though we lack a general procedure for establishing regularity results of this kind. The situation with supercritical equations is far more subtle. To start with, an equation that we now call supercritical22 may in fact turn out to be critical, or even subcritical, upon the discovery of additional a priori estimates. Thus an important question concerning the issue of criticality, and consequently that of singular behavior, is: are there other, stronger, local a priori bounds that cannot be derived from Noether’s principle? The discovery of such a bound would be a major event in both mathematics and physics.

Once we understand that the presence of singularities in our basic evolution equations is unavoidable, we have to face the question of whether they can somehow be accommodated by a more general concept of what a solution is or whether their structure is such that the equation itself, indeed the physical theory that it underlies, becomes meaningless. An acceptable concept of a generalized solution should, of course, preserve the deterministic nature of the equations: in other words, it should be uniquely determined from its Cauchy data.

Finally, once an acceptable concept of generalized solutions is found, we would like to use it to determine some important qualitative features, such as longterm asymptotic behavior. One can formulate a limitless number of such questions, the answers to which will vary from equation to equation.

(ii) Understand in a rigorous mathematical fashion the range of validity of various approximations. The equations obtained by various limiting procedures or phenomenological assumptions can of course be studied in their own right, as the examples that we have referred to above are. However, they present us with additional problems to do with the mechanics of how they are derived from equations that we regard as more fundamental. It is entirely possible, for example, that the dynamics of a derived system of equations leads to behavior that is incompatible with the assumptions made in its derivation. Alternatively, a particular simplifying assumption, such as spherical symmetry in general relativity or zero vorticity for compressible fluids, may turn out to be unstable at large scales and therefore not a reliable predictor of the general case. These and other similar situations lead to important dilemmas: should we persist in studying the approximate equations even when, in many cases, we face formidable mathematical difficulties (some which may turn out to be quite pathological and are perhaps related to the nature of the approximation), or should we abandon them in favor of the original system or a more suitable approximation? Whatever one may feel about this in any specific situation, it is clear that the problem of understanding, rigorously, the range of validity of various approximations is one of the fundamental goals in PDEs.

(iii) Devise and analyze the right equation for studying the specific geometric or physical problem at hand. This last goal is equally important even though it is necessarily vague. The enormously important role played by PDEs in various branches of mathematics is more evident than ever. One looks in awe at how equations such as the Laplace, heat, wave, Dirac, KdV, Maxwell, Yang– Mills, and Einstein equations, which were originally introduced in specific physical contexts, turned out to have very deep applications to seemingly unrelated problems in areas such as geometry, topology, algebra, and combinatorics. Other PDEs appear naturally in geometry when we look for embedded objects with optimal geometric shapes, such as solutions to isoperimetric problems, minimal surfaces, surfaces of least distortion or minimal curvature, or, more abstractly, connections, maps, or metrics with distinguished properties. They are variational in character, just like the main equations of mathematical physics. Other equations have been introduced with the goal of allowing one to deform a general object, such as a map, connection, or metric, to an optimal one. They usually arise in the form of geometric, parabolic flows. The most famous example of this is Ricci flow, first introduced by Richard Hamilton, who hoped to use it to deform Riemannian metrics into Einstein metrics. Similar ideas were used earlier to construct, for example, stationary harmonic maps with the help of a harmonic heat flow, and self-dual Yang–Mills connections with the help of a Yang–Mills flow. In addition to the successful use of Ricci flow to settle the Poincaré conjecture in three dimensions, another remarkable recent example of the usefulness of geometric flows is that of the inverse mean flow, first introduced by Geroch, to settle the so-called Riemannian version of the Penrose inequality.

# Further Reading

Brezis, H., and F. Browder. 1998. Partial differential equations in the 20th century. Advances in Mathematics 135: 76–144.   
Constantin, P. 2007. On the Euler equations of incompressible fluids. Bulletin of the American Mathematical Society 44:603–21.   
Evans, L. C. 1998. Partial Differential Equations. Graduate Studies in Mathematics, volume 19. Providence, RI: American Mathematical Society.   
John, F. 1991. Partial Differential Equations. New York: Springer.   
Klainerman, S. 2000. PDE as a unified subject. In \*GAFA 2000\*, Visions in Mathematics—Towards 2000 (special issue of Geometric and Functional Analysis), part 1, pp. 279–315.   
Wald, R. M. 1984. General Relativity. Chicago, IL: Chicago University Press.

# IV.13 General Relativity and the Einstein Equations

Mihalis Dafermos

Einstein’s formulation of general relativity represents one of the great triumphs of modern physics and provides the currently accepted classical theory that unifies gravitation, inertia, and geometry. The Einstein equations are the mathematical embodiment of this theory.

The definitive form of the equations,

$$
R _ {\mu \nu} - \frac {1}{2} R g _ {\mu \nu} = 8 \pi T _ {\mu \nu}, \tag {1}
$$

was attained in November 1915; this was the final act of Einstein’s eight-year struggle to generalize his principle of relativity so as to encompass gravitation, which had been described in the earlier “Newtonian” theory by the Poisson equation

$$
\frac {\partial^ {2} \phi}{\partial x ^ {2}} + \frac {\partial^ {2} \phi}{\partial y ^ {2}} + \frac {\partial^ {2} \phi}{\partial z ^ {2}} = 4 \pi \mu \tag {2}
$$

for the potential φ and mass density μ.

An obvious contrast between the Einstein equations (1) and the Poisson equation (2) is that the mysterious notation of the former makes it far less obvious what they even mean. This has given the subject of general relativity a reputation for difficulty and impenetrability. However, this reputation is to some extent unwarranted. Both (1) and (2) represent the culmination of revolutionary theories whose formulations presuppose a complicated conceptual framework. For better or for worse, however, the structure necessary to formulate Poisson’s equation has been incorporated into our traditional mathematical notation and school education. As a result, R3, with its Cartesian coordinate system, and notions such as functions, partial derivatives, masses, forces, and so on, are familiar to people with a general mathematical background, while the conceptual structure of general relativity is much less so, both with respect to its basic physical notions and with respect to the mathematical objects that are needed to model them. However, once one comes to terms with these, the equations turn out to be more natural and, one might even dare say, simpler.

Thus, the first task of this article is to explain in more detail the conceptual structure of general relativity. Our aim will be to make it clear what the equations (1) actually denote, and, moreover, why they are in a certain sense the simplest equations one can write down, given the general framework of the theory. This in turn will require us to review special relativity and its implications for the structure of matter, which will bring us to the unified concept of stress–energy–momentum, described by a tensorial object T . Finally, we will join Einstein in his inspired leap to the notion of a general four-dimensional Lorentzian manifold (M, g) that represents our space-time continuum. We shall see that equation (1) expresses a relationship between the tensor T and the geometry of g as expressed in its so-called curvature.

There is more to truly understanding a theory than merely knowing how to write down its governing equations. General relativity is associated with some of the most spectacular predictions of twentieth-century physics: gravitational collapse, black holes, space-time singularities, the expansion of the universe. These phenomena (which were completely unknown in 1915 and thus played no role in the formulation of the equations (1)) revealed themselves only when the conceptual issues surrounding the problem of global dynamics of solutions were understood. This took a surprisingly long time, though the story is not as well-known as the heroic struggle to attain (1). The article will conclude with a very brief glimpse into the fascinating dynamics of the Einstein equations.

# 1 Special Relativity

# 1.1 Einstein, 1905

Einstein’s 1905 formulation of special relativity stipulated that all fundamental laws of physics should be invariant under Lorentz transformations of the frame of reference defined by x, y, z, and t. A Lorentz transformation is any composition of translations, rotations, and the Lorentz boost, which is given by the formulas

$$
\left. \begin{array}{l} \tilde {x} = \frac {x - v t}{\sqrt {1 - v ^ {2} / c ^ {2}}}, \quad \tilde {y} = y, \\ \tilde {t} = \frac {t - v x / c ^ {2}}{\sqrt {1 - v ^ {2} / c ^ {2}}}, \quad \tilde {z} = z, \end{array} \right\} \tag {3}
$$

where c is a certain constant and $| \nu | < c .$ Thus, Einstein’s stipulation was that if one changes coordinates by means of a Lorentz transformation, then the form of all fundamental equations will remain the same. This set of transformations had already been identified in the context of the study of the vacuum Maxwell equations for the electric field E and magnetic field B:

$$
\left. \begin{array}{c c} \nabla \cdot \boldsymbol {E} = 0, & \nabla \cdot \boldsymbol {B} = 0, \\ c ^ {- 1} \partial_ {t} \boldsymbol {B} + \boldsymbol {\nabla} \times \boldsymbol {E} = 0, & c ^ {- 1} \partial_ {t} \boldsymbol {E} - \boldsymbol {\nabla} \times \boldsymbol {B} = 0. \end{array} \right\} \tag {4}
$$

Indeed, the Lorentz transformations are precisely the transformations that keep the form of the above equations invariant if we also transform E and B appropriately. Their significance was emphasized by poincaré [VI.61]. However, it was Einstein’s profound insight to elevate this invariance to the status of fundamental physical principle, despite its incompatibility with what we now usually call Galilean relativity, which corresponds to taking $c \quad  \quad \infty$ in (3). A surprising consequence of Lorentz invariance is that the notion of simultaneity is not absolute but depends on the observer: given two distinct events that occur at $( t , x , y , z )$ and $( t , x ^ { \prime } , y ^ { \prime } , z ^ { \prime } )$ , it is easy to find a Lorentz transformation such that the transformed events no longer have the same t-coordinate.

It follows from a celebrated result in partial differential equations known as the strong Huygens principle, applied to (4), that electromagnetic disturbances in vacuum propagate with speed $^ { c , }$ which we thus identify as the speed of light. In view of Lorentz invariance, this statement is independent of the frame! A further postulate of the principle of relativity is that physical theories should not allow massive particles to move at speeds (as measured in any frame) greater than or equal to c.

# 1.2 Minkowski, 1908

Einstein’s understanding of special relativity was “algebraic.” It was minkowski [VI.64] who first understood its underlying geometric structure, namely, that the content of the principle was contained in the metric element

$$
- c ^ {2} \mathrm{d} t ^ {2} + \mathrm{d} x ^ {2} + \mathrm{d} y ^ {2} + \mathrm{d} z ^ {2} \tag {5}
$$

defined on $\mathbb { R } ^ { 4 }$ with coordinates $( t , x , y , z )$ . We call $\mathbb { R } ^ { 4 }$ endowed with the metric (5) Minkowski space-time and denote it $\mathbb { R } ^ { 3 + 1 }$ . Points of $\mathbb { R } ^ { 3 + 1 }$ are referred to as events. The expression (5) is classical notation for the inner product defined on tangent vectors $\nu \_ =$ $( c ^ { - 1 } \nu ^ { 0 } , \nu ^ { 1 } , \nu ^ { 2 } , \nu ^ { 3 } )$ , $\textbf { \em w } = \mathbf { \phi } ( c ^ { - 1 } w ^ { 0 } , w ^ { 1 } , w ^ { 2 } , w ^ { 3 } )$ on $\mathbb { R } ^ { 4 }$ by

$$
\langle \boldsymbol {v}, \boldsymbol {w} \rangle = - v ^ {0} w ^ {0} + v ^ {1} w ^ {1} + v ^ {2} w ^ {2} + v ^ {3} w ^ {3}. \tag {6}
$$

The Lorentz transformations constitute precisely the symmetry group of the geometry defined by (5). Einstein’s principle of relativity could now be understood as the principle that the fundamental equations of physics must refer to space-time only through geometric quantities: that is, quantities that can be defined purely in terms of the metric. For example, from this point of view the reason that the notion of absolute simultaneity is not allowed is that it depends on a privileged hyperplane through any given point of $\mathbb { R } ^ { 3 + 1 }$ . But there are Lorentz transformations that preserve the metric and send this hyperplane to another one through the given point, so nothing in the metric can pick out one particular hyperplane. Note that if a physical theory makes use of geometric quantities only, then it is automatically invariant under Lorentz transformations: this observation renders many complicated calculations unnecessary.

Let us explore this geometric point of view further. Note that nonzero vectors v are naturally classified by the inner product $\langle \cdot , \cdot \rangle$ into three types, called timelike, null, and spacelike, according to whether $\langle \pmb { \nu } , \pmb { \nu } \rangle <$ $0 , ~ \langle \pmb { \nu } , \pmb { \nu } \rangle ~ = ~ 0 , ~ \mathrm { o r } ~ \langle \pmb { \nu } , \pmb { \nu } \rangle ~ > ~ 0$ , respectively. Idealized point particles traverse curves γ through space-time; these are called the world lines of the corresponding particles. The postulate (referred to earlier) that speed in any frame of reference is bounded by the speed of light c can now be formulated as the following statement: if γ is the world line of a particle, then the vector $\mathrm { d } \pmb { \gamma } / \mathrm { d } s$ must be timelike. (Null lines correspond to light rays in the geometric optics limit of (4).) This statement is independent of the parameter s of $\mathbf { \delta } _ { \mathbf { \mathcal { X } } , \mathbf { \delta } }$ but for world lines we shall always assume that dt/ds $> 0 .$ . To phrase this more geometrically, $\langle \mathrm { d } { \pmb { y } } / \mathrm { d } s , ( c ^ { - 1 } , 0 , 0 , 0 ) \rangle \ < \ 0$ , which we interpret as the statement that γ is future-directed.

We can now define the “length” of the world line of a particle by

$$
\begin{array}{l} L (\boldsymbol {\mathcal {Y}}) = \int_ {s _ {1}} ^ {s _ {2}} \sqrt {- \langle \dot {\boldsymbol {\mathcal {Y}}} , \dot {\boldsymbol {\mathcal {Y}}} \rangle} d s \\ = \int_ {s _ {1}} ^ {s _ {2}} \sqrt {c ^ {2} \left(\frac {\mathrm{d} t}{\mathrm{d} s}\right) ^ {2} - \left(\frac {\mathrm{d} x}{\mathrm{d} s}\right) ^ {2} - \left(\frac {\mathrm{d} y}{\mathrm{d} s}\right) ^ {2} - \left(\frac {\mathrm{d} z}{\mathrm{d} s}\right) ^ {2}} \mathrm{d} s. \tag {7} \\ \end{array}
$$

Classically, the above expression would have been written simply as

$$
L (\boldsymbol {\mathcal {Y}}) = \int_ {\boldsymbol {\mathcal {Y}}} \sqrt {- (- c ^ {2} \mathrm{d} t ^ {2} + \mathrm{d} x ^ {2} + \mathrm{d} y ^ {2} + \mathrm{d} z ^ {2})},
$$

which explains the notation (5). We refer to the quantity $c ^ { - 1 } L ( \pmb { \chi } )$ as proper time. This is the time that is relevant in local physical processes; in particular, if you are the particle traversing the world line γ, then $c ^ { - 1 } L ( \pmb { \chi } )$ is the time that you will feel.

The metric (5) contains three-dimensional Euclidean geometry

$$
\mathrm{d} x ^ {2} + \mathrm{d} y ^ {2} + \mathrm{d} z ^ {2},
$$

restricted to $t \ = \ 0 ,$ say. More interestingly, it also contains non-Euclidean geometry

$$
\left(1 - \frac {x}{r}\right) d x ^ {2} + \left(1 - \frac {y}{r}\right) d y ^ {2} + \left(1 - \frac {z}{r}\right) d z ^ {2}
$$

when it is restricted to the hypersurface $t = c ^ { - 1 } r =$ $c ^ { - 1 } \sqrt { x ^ { 2 } + y ^ { 2 } + z ^ { 2 } }$ . It is hard to overestimate how revolutionary the notion was that the time of physical processes (including our very sensations) and the length of measuring rods are two interdependent aspects of a geometric structure that naturally lives on a four-dimensional space-time continuum. Indeed, even Einstein initially rejected Minkowski space-time, preferring to retain the independent reality of a definite “space,” albeit a space with a relative notion of simultaneity. Only as a result of his search for general relativity did he realize that this view is fundamentally untenable. We shall return to this in section 3.

# 2 Relativistic Dynamics and the Unification of Energy, Momentum, and Stress

Besides the space-time concept and its geometrization, the principle of relativity led to a profound rearrangement and unification of the fundamental concepts of dynamics: mass, energy, and momentum. Einstein’s celebrated relation between mass and energy in the rest frame,

$$
E _ {0} = m c ^ {2}, \tag {8}
$$

is the best-known expression of one aspect of this unification. This relation arises naturally when one attempts to generalize Newton’s second law m(dv/dt)  f to a relation between 4-vectors in Minkowski space.

General relativity has to be formulated in terms of fields rather than particles. As a first step toward understanding it, let us look at continuous media. Now, instead of particles we consider matter fields; the unification of dynamical concepts encompasses what is known as stress, and its complete expression is embodied by the so-called stress–energy–momentum tensor T . This tensor is fundamental to general relativity, so we have no choice but to familiarize ourselves with it. It will be the key to the form of the Einstein equations (1) as well as to the object on their right-hand side.

For each point $\pmb q \in \mathbb { R } ^ { 3 + 1 }$ , the stress–energy–momentum tensor field T gives us a map

$$
\boldsymbol {T}: \mathbb {R} _ {\boldsymbol {q}} ^ {4} \times \mathbb {R} _ {\boldsymbol {q}} ^ {4} \rightarrow \mathbb {R} \tag {9}
$$

defined by the formula

$$
\boldsymbol {T} (\boldsymbol {w}, \tilde {\boldsymbol {w}}) = \sum_ {\alpha , \beta = 0} ^ {3} T _ {\alpha \beta} w ^ {\alpha} \tilde {w} ^ {\beta}.
$$

Here, $T _ { \alpha \beta } = T _ { \beta \alpha }$ for each α and β. By $\mathbb { R } _ { q } ^ { 4 }$ we mean the space of vectors at $\mathbf { \delta } _ { \mathbf { \delta } _ { q } } \mathbf { . }$ (In Minkowski coordinates, we often identify $\mathbb { R } ^ { 4 }$ with $\mathbb { R } _ { q } ^ { 4 }$ , but it will be important to distinguish between the two when considering arbitrary coordinates in section 3.2.) Bilinear maps of the form (9) are known as covariant 2-tensors.

If the only matter present is described by what is known as a perfect fluid, then the components of T are given by

$$
T _ {0 0} = (\rho + p) u ^ {0} u ^ {0} - p, \quad T _ {0 i} = (\rho + p) u ^ {i} u ^ {0},
$$

$$
T _ {i j} = (\rho + p) u ^ {i} u ^ {j} + p \delta^ {i j},
$$

where u is the 4-velocity, a timelike vector normalized such that $\langle { \pmb u } , { \pmb u } \rangle = - c ^ { 2 } , \rho$ is the mass–energy, p is the pressure, and where $\delta _ { i j } = 1 \mathrm { i f } i = j , 0 \mathrm { i f } i \neq j ,$ , and i and j range over 1, 2, 3. Greek indices will range over 0, 1, 2, 3. We identify T00 with energy, $T _ { 0 i }$ with momentum, and $T _ { i j }$ with stress. These notions are clearly framedependent. Finally, observe that $T ( \pmb { u } , \pmb { u } ) = \rho c ^ { 2 }$ . This is the field-theoretic version of the famous equation (8).

In general, T is derived from the totality of all the matter fields by constitutive functions that depend on the nature of the matter fields and their interactions. We need not worry here about such things. But, regardless of the nature of the matter fields involved, we always postulate that the following equations are satisfied:

$$
- \partial_ {0} T _ {0 \alpha} + \sum_ {i = 1} ^ {3} \partial_ {i} T _ {i \alpha} = 0.
$$

Defining $\nabla ^ { 0 } = - \partial _ { 0 } , \nabla ^ { i } = \partial _ { i }$ , and introducing the Einstein summation convention, under which summation is implicit when an index appears both upstairs and downstairs, we may rewrite this as

$$
\nabla^ {\mu} T _ {\mu \nu} = 0. \tag {10}
$$

These equations are Lorentz invariant.

The above relations embody the conservation $o f$ stress–energy–momentum at a differential level. Integrating (10) between homologous hypersurfaces and applying the Minkowski-space version of the divergence theorem, one obtains global balance laws. If one assumes that $T _ { \alpha \beta }$ is compactly supported, then, integrating between $t = t _ { 1 }$ and $t = t _ { 2 } ,$ , one obtains

$$
\int_ {t = t _ {2}} T _ {0 \alpha} \mathrm{d} x ^ {1} \mathrm{d} x ^ {2} \mathrm{d} x ^ {3} = \int_ {t = t _ {1}} T _ {0 \alpha} \mathrm{d} x ^ {1} \mathrm{d} x ^ {2} \mathrm{d} x ^ {3}. \tag {11}
$$

With respect to the chosen Lorentz frame, the zeroth component of the above equation represents the conservation of total energy, while the remaining components represent conservation of total momentum.

In the case of a perfect fluid, if we close the system (10) by adjoining a conservation law for particle number

$$
\nabla^ {\alpha} (n \pmb {u} _ {\alpha}) = 0
$$

and postulate constitutive relations between $\rho , p , { \mathrm { p a r } } -$ ticle number density n, and entropy per particle s, compatible with the laws of thermodynamics, then we arrive at the so-called relativistic Euler equations.

# 3 From Special to General Relativity

With the elements of special relativity at hand, together with their deep implications for the nature of energy, momentum, and stress, we can now pass to the formulation of general relativity.

# 3.1 The Equivalence Principle

Einstein understood as early as 1907 that the most profound aspect of the gravitational force could not be described within the relativity principle as he had formulated it in 1905. This aspect is what he called the equivalence principle.

The easiest setting in which to understand this principle is that of the “test particle” with velocity v(t) in a fixed gravitational field φ. In this case, we have that the classical gravitational force is given by $f = - m \nabla \phi$ , and we may rewrite Newton’s second law m(dv/dt) f as

$$
\frac {\mathrm{d} \boldsymbol {v}}{\mathrm{d} t} = - \nabla \phi . \tag {12}
$$

Notice that the mass m has dropped out! Thus, the gravitational field accelerates all objects at a given position in the same way. This explains the fact, recorded already in late antiquity by Ioannes Philoponus and popularized in Western Europe by Galileo, that the time it takes objects to fall from a given height is independent of their weight.

It was Einstein who first interpreted this property as a sort of covariance with respect to transformations to noninertial, that is to say accelerated, frames. For instance, in the case of a constant gravitational field, which corresponds to the case $\phi ( z ) = f z$ , we can pass to the accelerated frame

$$
\tilde {z} = z + \frac {1}{2} f t ^ {2}
$$

and write (12) as

$$
\frac {\mathrm{d} \boldsymbol {v}}{\mathrm{d} t} = 0. \tag {13}
$$

Similarly, one can reverse the argument to “simulate” a gravitational field when none is present by expressing (13) in an accelerated frame.

# 3.2 Vectors, Tensors, and Equations in General Coordinates

Exactly what the equivalence principle means in general is somewhat obscure and has been the subject of debate ever since Einstein introduced it. Nevertheless, the above considerations suggest that, even in the absence of gravity, it would be useful to know how various objects and equations appear when expressed in arbitrary coordinate systems. That is to say, let us change from our Minkowski coordinates $x ^ { 0 } , x ^ { 1 } , x ^ { 2 } , x ^ { 3 }$ to the most general coordinate system, which we shall write as $\bar { x } ^ { \bar { \mu } } = \bar { x } ^ { \bar { \mu } } ( x ^ { 0 } , x ^ { 1 } , x ^ { 2 } , x ^ { 3 } )$ , where $\bar { \mu }$ ranges over $0 , 1 , 2 , 3$ .

Expressing scalar functions in arbitrary coordinates poses no problem. But what about vector fields? If v is a vector field expressed in Minkowski coordinates as $( \nu ^ { 0 } , \nu ^ { 1 } , \nu ^ { 2 } , \nu ^ { 3 } )$ , how do we express v in our new coordinates $\bar { x } ^ { \bar { \mu } \gamma }$

One has to think a bit about what a vector field actually is. The correct point of view is to consider a vector field v as a first-order differential operator defined (using Einstein’s summation convention) by $\pmb { \nu } ( f ) \ =$ $\nu ^ { \mu } \partial _ { \mu } f$ . So we seek $\nu ^ { \bar { \mu } }$ such that $\pmb { \nu } ( f ) = \nu ^ { \bar { \mu } } \partial _ { \bar { \mu } } f$ for all functions $f .$ . The chain rule then gives us our answer:

$$
\nu^ {\bar {\mu}} = \frac {\partial \bar {x} ^ {\bar {\mu}}}{\partial x ^ {\nu}} \nu^ {\nu}. \tag {14}
$$

What about tensors, such as the stress–energy–momentum tensor ${ \pmb T } \smash { \ k _ { 0 } }$ In view of the definition (9), we seek

$T _ { \bar { \mu } \bar { \nu } }$ such that

$$
\boldsymbol {T} (\boldsymbol {u}, \boldsymbol {v}) = T _ {\bar {\mu} \bar {\nu}} u ^ {\bar {\mu}} v ^ {\bar {\nu}}, \tag {15}
$$

where the numbers $u ^ { \bar { \mu } }$ are the components of u with respect to the coordinates $\bar { x } ^ { \bar { \mu } }$ as we have just calculated them above. (Note that these components depend on the point q. This is why it is now essential to distinguish $\mathbb { R } _ { q } ^ { 4 }$ from $\mathbb { R } ^ { 4 } . )$ Again, the chain rule gives us the answer:

$$
T _ {\bar {\mu} \bar {\nu}} = T _ {\mu \nu} \frac {\partial x ^ {\nu}}{\partial \bar {x} ^ {\bar {\nu}}} \frac {\partial x ^ {\mu}}{\partial \bar {x} ^ {\bar {\mu}}}.
$$

Classically, we write

$$
\boldsymbol {T} = T _ {\bar {\mu} \bar {\nu}} \mathrm{d} \bar {x} ^ {\bar {\mu}} \mathrm{d} \bar {x} ^ {\bar {\nu}} = T _ {\mu \nu} \mathrm{d} x ^ {\mu} \mathrm{d} x ^ {\nu}.
$$

One can interpret the above as a shorthand notation for (15), but it also tells us how to compute $T _ { \bar { \mu } \bar { \nu } }$ from $T _ { \mu \nu }$ by formally applying the chain rule to ${ \mathrm { d } } \bar { x } ^ { \bar { \mu } }$ .

There is another covariant symmetric 2-tensor besides T that is relevant here. This is the Minkowski metric itself. Indeed, the classical form of the Minkowski metric (5) corresponds to the representation

$$
\eta_ {\mu \nu} \mathrm{d} x ^ {\mu} \mathrm{d} x ^ {\nu},
$$

where the $\eta _ { \mu \nu }$ for Minkowski coordinates $x ^ { \mu }$ are given by $\eta _ { 0 0 } = - 1 , \eta _ { 0 i } = 0 , \eta _ { i j } = 1 \mathrm { ~ i f ~ } i = j$ , and $\eta _ { i j } = 0$ i f $i \neq j .$ To avoid the cumbersome notation $\langle \cdot , \cdot \rangle$ , let us refer to the Minkowski metric as η. Following the above, we may express η in general coordinates $\bar { x } ^ { \bar { \mu } }$ by

$$
\eta_ {\bar {\mu} \bar {\nu}} \mathrm{d} \bar {x} ^ {\bar {\mu}} \mathrm{d} \bar {x} ^ {\bar {\nu}},
$$

where $\eta _ { \bar { \mu } \bar { \nu } }$ is computed by formal application of the chain rule.

It is clear that if one tries to transform an equation such as (10) into general coordinates, then the components of η and their derivatives will appear in the equations. Einstein (always thinking “algebraically”) was seeking laws of motion for both matter and the gravitational field that would have the same form in all coordinate systems. As he understood it, this meant that all objects that appear should transform as tensors and should be considered a priori “unknown.” He referred to this principle as “general covariance.” This suggests that η should be replaced by an unknown symmetric 2-tensor. Let us call this 2-tensor $\pmb { g } .$ One can of course try to write down an equation for the “unknown” $\pmb { g }$ that forces it to be the “known” Minkowski metric η. Thus, “general covariance” per se does not force one to abandon η. But in view of the fact that g and T have the same number of components, it was a natural step to consider g as the embodiment of the gravitational field and to try to look for an equation that related g and T directly. In this way, the framework of general relativity was born.

# 3.3 Lorentzian Geometry

The profound insight of replacing the fixed Minkowski η with a dynamic g brought Einstein to what we now call Lorentzian geometry. Lorentzian geometry generalizes Minkowski geometry following the blueprint of riemann [VI.49]. That is, we replace the Minkowski metric η by a general map

$$
\boldsymbol {g}: \mathbb {R} _ {\boldsymbol {q}} ^ {4} \times \mathbb {R} _ {\boldsymbol {q}} ^ {4} \to \mathbb {R}.
$$

In other words, we replace η by a symmetric covariant 2-tensor, which is expressed in arbitrary coordinates $x ^ { \mu }$ by

$$
g _ {\mu \nu} \mathrm{d} x ^ {\mu} \mathrm{d} x ^ {\nu}.
$$

Moreover, we require that at each point q the bilinear form $\mathbf { \sigma } _ { \mathbf { \sigma } _ { g } ( \cdot , \cdot ) }$ can be diagonalized to the Minkowski form (6). Loosely speaking, a Lorentzian metric is one that “looks locally like the Minkowski metric,” just as a riemannian metric [I.3 §6.10] looks locally like the Euclidean metric.

Just as with the Minkowski metric, the bilinear form g permits us to classify nonzero vectors $\nu _ { q }$ at a point q as timelike, null, or spacelike and to define proper times of world lines $\gamma ( s ) ~ = ~ ( x ^ { 0 } ( s ) , x ^ { 1 } ( s ) , x ^ { 2 } ( s ) , x ^ { 3 } ( s ) )$ by the formula $( 7 ) ,$ but with $\langle \dot { \pmb { y } } , \dot { \pmb { y } } \rangle$ replaced by $g _ { \mu \nu } \dot { x } ^ { \mu } \dot { x } ^ { \nu }$ . It is in this sense that we can speak of the geometry of ${ \pmb g } .$

In view of Minkowski’s formulation of the special relativity principle as the statement that the equations of physics refer to space-time only through geometric quantities associated with the Minkowski metric, it is natural to look for a generalization of this principle, and indeed a suitable version immediately suggests itself. It is the principle that the equations of physics refer to the space-time coordinates only via geometric quantities naturally associated with ${ \pmb g } .$ .

We saw earlier that the kinematic constraint on “test particles,” as formulated geometrically for the Minkowski metric, was that $\mathrm { d } \pmb { \gamma } / \mathrm { d } s$ should be timelike; this makes sense for an arbitrary Lorentzian metric. But how does one formulate differential equations? For instance, how does one formulate an analogue of (10) that refers only to $\pmb { g } \smash { ? }$

It turned out that in the Riemannian case, a set of natural geometric concepts suitable for the task had already been developed in the nineteenth and early twentieth centuries by Riemann, Bianchi, Christoffel, Ricci, and Levi-Civita. These carry over directly to the Lorentzian case.

One begins by defining the so-called Christoffel symbols $T _ { \mu \nu } ^ { \lambda }$ by

$$
\Gamma_ {\mu \nu} ^ {\lambda} = \frac {1}{2} g ^ {\lambda \rho} (\partial_ {\mu} g _ {\rho \nu} + \partial_ {\nu} g _ {\mu \rho} - \partial_ {\rho} g _ {\mu \nu}).
$$

Here, the numbers $g ^ { \mu \nu }$ are the components of the “inverse metric” of g: that is, they are the unique solution to the equation $g ^ { \mu \nu } g _ { \nu \lambda } ~ = ~ \delta _ { \lambda } ^ { \mu }$ where, as usual, $\delta _ { \lambda } ^ { \mu } = 1 \mathrm { i f } \lambda = \mu$ and 0 otherwise. (It turns out that $g ^ { \mu \nu }$ is very useful for the calculational gymnastics that are typical of tensor analysis when it exploits the Einstein summation convention.)

One can then define a differential operator $\nabla _ { \mu }$ called a connection, which acts on vector fields by

$$
\nabla_ {\mu} \boldsymbol {v} ^ {\nu} = \partial_ {\mu} \boldsymbol {v} ^ {\nu} + \Gamma_ {\mu \lambda} ^ {\nu} \boldsymbol {v} ^ {\lambda} \tag {16}
$$

and on covariant 2-tensors by

$$
\nabla_ {\lambda} T _ {\mu \nu} = \partial_ {\lambda} T _ {\mu \nu} - \Gamma_ {\lambda \mu} ^ {\sigma} T _ {\sigma \nu} - \Gamma_ {\lambda \nu} ^ {\sigma} T _ {\mu \sigma}. \tag {17}
$$

The left-hand sides of (16) and (17) define tensors that can be expressed in any coordinate system by a formal application of the chain rule.

With the help of this differential operator, one could now write the analogue of equations (10) for an arbitrary metric g as

$$
\nabla^ {\mu} T _ {\mu \nu} = 0, \tag {18}
$$

where $\nabla ^ { \mu } = g ^ { \mu \nu } \nabla \cdot$ ν refers to the connection associated with ${ \pmb g } .$ .

If we consider a limit as the matter field becomes concentrated at a point, or rather as the stress–energy– momentum tensor $T _ { \mu \nu }$ is nonzero only on a world line, then this curve will be a geodesic of $\pmb { g } \colon$ that is, a curve that locally maximizes the proper time defined by g. These are the analogues of straight timelike lines in Minkowski space. In this limit, the motion of the matter does not depend on the nature of the stress– energy–momentum tensor, but only on the geometry of the metric that defines geodesics. Thus, all objects fall in the same way. These considerations give a concrete realization to the equivalence principle in general relativity.

Finally, it is important to remark that for a general metric $^ { g , }$ the identity (18) does not imply global conservation laws (11) for “total energy” and “total momentum.” Such laws hold only if g has symmetries. The fact that the fundamental conservation laws survive in general only at the infinitesimal level is an important insight into the nature of these principles in physics.

# 3.4 Curvature and the Einstein Equations

It remains, then, to give a set of equations for the metric g that relate it to T . In anticipation of a Newtonian limit, we expect these equations to be second order, and we expect them to implement “general covariance” in the simplest way possible: they should refer to no other structure but g itself and T .

Again, Riemannian geometry provides ready-made tensorial objects that are invariantly associated with ${ \pmb g } .$ One can define the Riemann curvature tensor

$$
R _ {\mu \nu \lambda \rho} \mathrm{d} x ^ {\mu} \mathrm{d} x ^ {\nu} \mathrm{d} x ^ {\lambda} \mathrm{d} x ^ {\rho}
$$

with components given by

$$
R _ {\mu \nu \lambda \rho} = g _ {\mu \sigma} (\partial_ {\rho} \Gamma_ {\nu \lambda} ^ {\sigma} - \partial_ {\lambda} \Gamma_ {\nu \rho} ^ {\sigma} + \Gamma_ {\nu \lambda} ^ {\tau} \Gamma_ {\tau \rho} ^ {\sigma} - \Gamma_ {\nu \rho} ^ {\tau} \Gamma_ {\tau \lambda} ^ {\sigma}).
$$

One can also define the Ricci curvature

$$
R _ {\mu \nu} \mathrm{d} x ^ {\mu} \mathrm{d} x ^ {\nu},
$$

a covariant symmetric 2-tensor with components given by

$$
R _ {\mu \nu} = g ^ {\lambda \rho} R _ {\mu \nu \lambda \rho},
$$

and the scalar curvature

$$
R = g ^ {\mu \nu} R _ {\mu \nu}.
$$

If g were the induced (Riemannian) metric on a 2- surface in $\mathbb { R } ^ { 3 }$ , then R would just be twice the Gauss curvature K. The above expressions should be thought of as complicated tensorial generalizations of Gauss curvature to several dimensions.

The final piece of the puzzle for the formulation of the Einstein equations (1) is provided by the following constraint that Einstein demanded: whatever the equation relating the metric and the stress–energy–momentum tensor of matter, (18) (the infinitesimal conservation of stress–energy–momentum) should hold as a consequence. Now, it turns out that for any metric $^ { g , }$ the so-called Bianchi identities imply that

$$
\nabla^ {\mu} (R _ {\mu \nu} - \frac {1}{2} g _ {\mu \nu} R) = 0. \tag {19}
$$

It is thus natural to postulate a linear relation between $T _ { \mu \nu }$ and the tensor $R _ { \mu \nu } - { \textstyle \frac { 1 } { 2 } } g _ { \mu \nu } R$ . The form

$$
R _ {\mu \nu} - \frac {1}{2} g _ {\mu \nu} R = 8 \pi G c ^ {- 4} T _ {\mu \nu} \tag {20}
$$

is then uniquely determined by the requirement that it should give the correct Newtonian limit when one makes the identifications

$$
g _ {0 0} \sim 1 + 2 \phi / c ^ {2}, \quad g _ {0 j} \sim 0, \quad g _ {i j} \sim (1 - 2 \phi / c ^ {2}) \delta_ {i j}.
$$

The form (1) corresponds to the usual units $G = c = 1$ . Note that (1), when written out explicitly, is nonlinear in the metric components $g _ { \mu \nu }$ .

Einstein did not stop at the Newtonian limit. By considering geodesic motion in solutions of the linearized equations (20), Einstein was able to determine the correct value for the anomalous precession of the perihelion of Mercury, an effect that Newtonian theory was unable to explain. Since (20) had no adjustable parameters after determining the Newtonian limit, this was a genuine test of the theory. A few years later the gravitational “bending” of light was observed. This had been calculated theoretically in the context of the geometric optics approximation where light rays follow null geodesics in a fixed space-time background. Post-Newtonian predictions of (1) have now been verified by various solar system tests, confirming general relativity in this regime to a high degree of accuracy.

One special case of (20) is when we postulate that $T _ { \mu \nu } = 0$ . The equations then simplify to

$$
R _ {\mu \nu} = 0. \tag {21}
$$

These are known as the vacuum equations. The Minkowski metric (5) is a particular solution (but not the only one!).

The vacuum equations can be derived formally as the euler–lagrange equations [III.94] corresponding to the so-called Hilbert Lagrangian:

$$
\mathcal {L} (\pmb {g}) = \int R \sqrt {- g} \mathrm{d} x ^ {0} \mathrm{d} x ^ {1} \mathrm{d} x ^ {2} \mathrm{d} x ^ {3}.
$$

(The expression ${ \sqrt { - g } } \mathrm { d } x ^ { 0 } \mathrm { d } x ^ { 1 } \mathrm { d } x ^ { 2 } \mathrm { d } x ^ { 3 }$ denotes the natural volume form associated with $\mathbf { _ { g . } } )$ hilbert [VI.63], who was following closely Einstein’s struggle to formulate a theory of gravity with a dynamic metric g, arrived at his Lagrangian (actually a more general version of the above yielding the coupled Einstein–Maxwell system) very shortly before Einstein obtained the general equations (20).

Many of the most interesting phenomena that come from the equations (20) are already present in the vacuum case (21). This is somewhat ironic, because it was the forms of T and (10) that dictated (20). Note, in contrast, that in the Newtonian theory (2), the “vacuum” equations $\mu = 0$ and standard boundary conditions at infinity imply $\phi = 0$ . Thus, the Newtonian theory of the vacuum is trivial.

The part of the curvature tensor $R _ { \mu \nu \lambda \rho }$ that is not forced to vanish from (21) is known as the Weyl curvature. This curvature measures the “tidal” distortion of families of geodesics. Thus, the “local strength” of gravitational fields in vacuum regions is related in the Newtonian limit to the tidal forces on macroscopic test matter, not the norm of the gravitational force.

# 3.5 The Manifold Concept

We have been able to get this far without really addressing the question of where the metric g is defined. In passing from the Minkowski metric to a general g, Einstein did not originally have in mind replacing the domain $\mathbb { R } ^ { 4 }$ . But it is clear in the Riemannian case from the theory of surfaces that the natural object for a metric to live on is not necessarily $\mathbb { R } ^ { 2 }$ but a general surface. For instance, the metric $\mathrm { d } \theta ^ { 2 } + \sin \theta \mathrm { d } \phi ^ { 2 }$ naturally lives on the sphere $\mathbb { S } ^ { 2 }$ . In saying this, we are to understand that one requires several coordinate systems of the type $( \theta , \phi )$ to cover all of $\mathbb { S } ^ { 2 }$ . The ndimensional generalization of the object where Riemannian or Lorentzian metrics naturally live is a manifold [I.3 §6.9]. Manifolds are the structures obtained by consistently smoothly pasting together local coordinate systems.

Thus, general relativity allows the space-time continuum not to be $\mathbb { R } ^ { 4 }$ but instead to be a general manifold ${ \mathcal { M } } ,$ which may very well be topologically inequivalent to $\mathbb { R } ^ { 4 }$ , just as $\mathbb { S } ^ { 2 }$ is inequivalent to $\mathbb { R } ^ { 2 }$ . We call the pair $( \mathcal { M } , \pmb { g } )$ a Lorentzian manifold. Properly put, the unknown in the Einstein equations is not just g but the pair $( \mathcal { M } , \pmb { g } )$ .

It is interesting that this fundamental fact, namely that the topology of space-time is not a priori determined by the equations, arises almost as an afterthought. Moreover, it was a thought that took many years to be clarified.

# 3.6 Waves, Gauges, and Hyperbolicity

When written out explicitly in arbitrary coordinates (try it!), the Einstein equations do not appear to be of any usual type, such as elliptic (like the poisson equation [IV.12 §1]), parabolic (like the heat equation [I.3 §5.4]), or hyperbolic (like the wave equation [I.3 §5.4]; see [IV.12 §2.5] for more about these different classes of PDEs). This is related to the fact that, given a solution, one can form a “new” solution by composing the old solution with a coordinate transformation. We can do this for new coordinate systems whose coordinate transformations differ from the identity only in a ball. This fact, known as the hole argument, confused Einstein and his mathematical collaborator Marcel Grossmann, who were thinking algebraically in terms of the form of the equations in coordinates, and temporarily led them to reject “general covariance.” The resulting backtracking delayed the final correct formulation of (1) by about two years. The geometric interpretation of the theory immediately suggests the resolution to the dilemma: such solutions are to be considered “the same” because they are the same from the point of view of all geometric measurements. In modern language, a solution to the Einstein vacuum equations (say) is an equivalence class [I.2 §2.3] of spacetimes (M, g), where two space-times are equivalent if there exists a diffeomorphism φ between them such that in any open set the metric has the same coordinate form when one identifies local coordinates by φ.

It turns out that once these conceptual issues are overcome, the Einstein equations can be viewed as hyperbolic. The easiest way to do this is to impose a gauge: that is to say, a certain restriction on the coordinate system. Specifically, one requires the coordinate functions $x ^ { \alpha }$ to satisfy the wave equation $\begin{array} { r } { \boxed { g ^ { } x ^ { \alpha } } = 0 , } \end{array}$ , where the d’Alembertian operator is defined by the formula

$$
\Box_ {g} = \frac {1}{\sqrt {- g}} \partial_ {\mu} (\sqrt {- g} g ^ {\mu \nu} \partial_ {\nu}).
$$

Such coordinates always exist locally and they are traditionally called harmonic coordinates, although the term wave coordinates would perhaps be more appropriate. The Einstein equation can then be written as a system

$$
\Box_ {g} g _ {\mu \nu} = N _ {\mu \nu} (\{g _ {\alpha \beta} \}, \{\partial_ {\gamma} g _ {\alpha \beta} \}),
$$

where $N _ { \mu \nu }$ is a nonlinear expression that is quadratic in the $\partial _ { \gamma } g _ { \alpha \beta }$ . In view of the Lorentzian signature of the metric, the above system constitutes what is known as a second-order nonlinear (but quasilinear) hyperbolic system.

At this point, it is instructive to make a comparison with the Maxwell equations. Suppose we are given an electric field E and a magnetic field B defined on Minkowski space. A 4-potential is a vector field A such that $E _ { i } = - \partial _ { i } A _ { 0 } - c ^ { - 1 } \partial _ { t } A _ { i }$ , and $\begin{array} { r } { B _ { i } = \sum _ { j , k = 1 } ^ { 3 } \epsilon _ { i j k } \partial _ { j } A _ { k } . } \end{array}$ (Here $\epsilon _ { 1 2 3 } = 1$ , and $\epsilon _ { i j k }$ is totally antisymmetric, i.e., it transforms to its negative under permutation of any two indices.) If one wishes to view A as the fundamental physical object, then one notices that if A is replaced by the field A˜, defined by the formula

$$
\tilde {\boldsymbol {A}} = \boldsymbol {A} + (- c ^ {- 1} \partial_ {t} \psi , \partial_ {1} \psi , \partial_ {2} \psi , \partial_ {3} \psi),
$$

where ψ is an arbitrary function, then A˜ is also a 4- potential for E and B. One can expect a determined equation for A only if one imposes further conditions on it: that is, if one “fixes the gauge.” (The terminology “gauge” is originally due to weyl [VI.80].) In the so-called Lorentz gauge

$$
\nabla^ {\mu} A _ {\mu} = 0,
$$

the Maxwell equations can be written

$$
\Box A _ {\mu} = - c ^ {- 2} \partial_ {t} ^ {2} A _ {\mu} + \sum_ {i} \partial_ {i} ^ {2} A _ {\mu} = 0,
$$

from which the wave properties are completely manifest. The gauge-symmetric point of view lived on to later twentieth century glory: the Yang–Mills equations, which are a nonlinear generalization of the Maxwell equations with a similar gauge symmetry, are the central part of the so-called standard model for particle physics.

The hyperbolicity property of the Einstein equations has two important repercussions. The first is that there should exist gravitational waves. This was noted by Einstein at least as early as 1918, essentially as a result of a linearized version of the considerations in the above discussion. The second is that there is a well-posed initial value problem [IV.12 §2.4] for the Einstein equations (1) with the domain-of-dependence property, when these are coupled with appropriate matter equations. In particular, this is true in the vacuum case (21). The proper conceptual framework to formulate the latter problem took a long time to get right, and was only completely understood through work of Choquet-Bruhat and Geroch in the 1950s and 1960s, based on the fundamental concept of global hyperbolicity due to Leray. Well-posedness means that one could associate a unique solution (in the vacuum case, a Lorentzian 4-manifold $( \mathcal { M } , g )$ satisfying (21)) with a suitable notion of initial data. Of course, “initial data” does not mean “data at time $t = 0 , "$ since the concept of t 0 is not geometric. Instead, the data take the form of some Riemannian 3-manifold $( \Sigma , \bar { g } )$ with a symmetric covariant 2-tensor K. The triple (Σ, g, K) ¯ has to satisfy the so-called Einstein constraint equations. But with this notion, the fundamental problem of general relativity, despite its revolutionary conceptual structure, is thoroughly classical: to determine the relation of the solution to initial data, that is to say, to determine the future from knowledge of the “present.” This is the problem of dynamics.

# 4 The Dynamics of General Relativity

In this final section we give a taste of our current mathematical understanding of the dynamics of the Einstein equations.

# 4.1 Stability of Minkowski Space and the Nonlinearity of Gravitational Radiation

In any physical theory in which one can formulate the problem of dynamics, the most basic question is the stability of the trivial solution. In other words, if we make a small change to the “initial conditions,” will the resulting change to the solution be small as well? In the case of general relativity, this is the question of stability of the Minkowski space-time $\mathbb { R } ^ { 3 + 1 }$ . This fundamental result was proven for the vacuum equations (21) in 1993 by Christodoulou and Klainerman.

The proof of the stability of Minkowski space made it possible to formulate the laws of gravitational radiation rigorously. Gravitational radiation is yet to be observed directly, but it has been inferred, originally by Hulse and Taylor, from the energy loss of a binary system. This work gave them the only Nobel prize (1993) directly associated with the Einstein equations! The blueprint for the mathematical formulation of the radiation problem is based on work of Bondi and later Penrose. One associates with the space-time $( \mathcal { M } , g )$ an ideal boundary “at infinity,” known as null infinity and denoted ${ \mathcal { I } } ^ { + }$ . Physically, the points of ${ \mathit { 1 } } ^ { + }$ correspond to observers who are far away from the isolated selfgravitating system but who are receiving its signals. Gravitational radiation can be identified with certain tensors defined on ${ \mathit { 1 } } ^ { + }$ from rescaled boundary limits of various geometric quantities. As Christodoulou was to discover, the laws of gravitational radiation are themselves nonlinear, and the nonlinearity is potentially relevant for observation.

# 4.2 Black Holes

Perhaps no prediction of general relativity is better known today than that of black holes.

The story of black holes begins with the so-called Schwarzschild metric:

$$
\begin{array}{l} - \left(1 - \frac {2 m}{r}\right) \mathrm{d} t ^ {2} + \left(1 - \frac {2 m}{r}\right) ^ {- 1} \mathrm{d} r ^ {2} \\ + r ^ {2} (\mathrm{d} \theta^ {2} + \sin^ {2} \theta \mathrm{d} \phi^ {2}). \tag {22} \\ \end{array}
$$

The parameter m here is a positive constant. This is a solution of the vacuum Einstein equations (21) that was found in 1916. The original interpretation of (22) was that it modeled the gravitational field in a vacuum region outside a star. That is to say, (22) was considered only in some coordinate range $r > R _ { 0 }$ , for an $R _ { 0 } > 2 m$ , and the metric was matched at $r = R _ { 0 }$ to a “static” interior metric satisfying the coupled Einstein–Euler system in the coordinate range $r \leqslant R _ { 0 }$ . (This latter metric is again of the form (22), but with $m = m ( r )$ such that $m  0 \mathrm { \ a s \ } r  0 . )$

From the theoretical point of view, a natural problem poses itself. Suppose we do away with the star altogether and try to consider (22) for all values of r . What happens then to the metric (22) at $r = 2 m ?$ In the (r , t) coordinates, the metric element appears to be singular. But this turns out to be an illusion! By a simple change of coordinates, one can easily extend the metric regularly as a solution of (21) beyond $r \ = \ 2 m$ . That is, there exists a manifold  that contains both a region $r > 2 m$ and a region $0 < r < 2 m$ , separated by a regular (null) hypersurface ${ \mathcal { H } } ^ { + }$ . The metric element (22) is valid everywhere except on ${ \mathcal { H } } ^ { + }$ , where it must be rewritten in regular coordinates.

It turns out that the hypersurface ${ \mathcal { H } } ^ { + }$ can be characterized by an exceptional global property: it defines the boundary of the region of space-time that can send signals to null infinity I+, or, in the physical interpretation, to distant observers. In general, the set of points that cannot send signals to null infinity ${ \mathit { 1 } } ^ { + }$ is known as the black hole region of space-time. Thus, the region $0 < r < 2 m$ is the black hole region of ${ \mathcal { M } } ,$ and ${ \mathcal { H } } ^ { + }$ is known as the event horizon.

These issues took a long time to be sorted out, partly because the language of global Lorentzian geometry was developed long after the original formulation of the Einstein equations. The global geometry of the extended space-time  was clarified by Synge in around 1950 and finally by Kruskal in 1960. The name “black hole” is due to the imaginative physicist John Wheeler. From their beginnings as a theoretical curiosity, black holes have become part of the accepted astrophysical explanation for a wide variety of phenomena, and in particular are thought to represent the end-state for the gravitational collapse of many stars.

# 4.3 Space-Time Singularities

A second natural problem poses itself in relation to the Schwarzschild metric (22), now considered in the region $r \mathrm { ~ < ~ }$ 2m of the extended space-time M: what happens at $r = 0 ?$

A computation reveals that as $r  0 ,$ , the Kretchmann scalar $R _ { \mu \nu \lambda \rho } R ^ { \mu \nu \lambda \rho }$ blows up. Since this expression is a geometric invariant, it follows that, unlike the situation at $r = 2 m$ , the space-time is not regularly extendable beyond 0. Moreover, timelike geodesics (freely falling observers in the test particle approximation) entering the black hole region reach $r = 0$ in finite proper time, so they are “incomplete” in the sense that they cannot be continued indefinitely. They thus “observe” the breakdown of the geometry of the space-time metric. Moreover, macroscopic observers approaching $r \ = \ 0$ are torn apart by the gravitational “tidal forces.”

In the early years of the subject, it was thought that this seemingly pathological behavior was connected to the high degree of symmetry of the Schwarzschild metric and that “generic” solutions would not exhibit such phenomena. That this is not the case was shown by Penrose’s celebrated incompleteness theorem of 1965. This states that solutions to the initial value problem for the Einstein equations coupled to appropriate matter will always contain such incomplete timelike or null geodesics if the initial data hypersurface is noncompact and contains what is known as a closed trapped surface. The Schwarzschild case may appear to suggest that such incomplete geodesics are associated with the curvature blowing up. However, the situation can in fact be very different, as is apparent in the celebrated Kerr solutions, a remarkable two-parameter family of solutions to the vacuum equations (21), discovered only in 1963, which are rotating versions of (22). In the Kerr solutions, incomplete timelike geodesics meet a so-called Cauchy horizon, a smooth boundary of the region of space-time that is uniquely determined by initial data.

The theorem of Penrose gives rise to two important conjectures. The first, known as weak cosmic censorship, says roughly that for generic physically plausible initial data for suitable Einstein-matter systems, geodesic incompleteness, if it occurs, is always confined to black hole regions. The second, strong cosmic censorship, says roughly that for generic admissible initial data, incompleteness of the solution is always associated with a local obstruction to extendability, such as the blow-up of curvature. The latter conjecture would ensure that the unique solution of the initial value problem is the only classical space-time that can arise from the data. That is to say, it would imply that classical determinism holds for the Einstein equations.

Both conjectures are false if we drop the assumption that the initial data are generic, and this is one reason for their difficulty. Indeed, Christodoulou has constructed spherically symmetric solutions of the coupled Einstein-scalar field system (arising from regular initial data) that are geodesically incomplete but do not contain black hole regions. Such space-times are said to contain naked singularities.

Naked singularities are easy to construct if one does not require that they arise from the collapse of regular initial data. An example is the Schwarzschild metric (22) for $m < 0 .$ . This metric, however, does not admit a complete asymptotically flat Cauchy hypersurface. This fact is related to the celebrated positive energy theorem of Schoen and Yau.

# 4.4 Cosmology

The space-times $( \mathcal { M } , \pmb { g } )$ discussed previously are all idealized representations of isolated systems. The “rest of the universe” is excised and replaced by an “asymptotically flat end”; far-away observers are placed at an ideal boundary “at infinity.” But what if we are more ambitious and consider our space-time $( \mathcal { M } , \pmb { g } )$ as representing the whole universe? The study of this latter problem is known as cosmology.

Observations suggest that on very large scales the universe is approximately homogeneous and isotropic. This is sometimes known as the Copernican principle. Interestingly, one cannot solve the Poisson equation (2) with a constant φ and constant nonzero μ on R4. Thus, in Newtonian physics, cosmology never became a rational science.1 General relativity, on the other hand, does admit homogeneous and isotropic solutions as well as their perturbations. Indeed, cosmological solutions of the Einstein equations were studied by Einstein himself, de Sitter, Friedmann, and Lemaitre in the early years of the subject.

When general relativity was formulated, the prevailing view was that the universe should be static. This led Einstein to add a term $\varLambda g _ { \mu \nu }$ to the left-hand side of his equations, fine-tuned so as to allow for such a solution. The constant Λ is known as the cosmological constant. The expansion of the universe is now considered to be an observational fact, beginning with the fundamental discoveries of Hubble. Expanding universes can be modeled to a first approximation by so-called Friedmann–Lemaitre solutions to the Einstein–Euler system, with various values of Λ. In the past direction, these solutions are singular: this singular behavior is often given the suggestive name “the big bang.”

# 4.5 Future Developments

The plethora of exact solutions of the Einstein equations gives us a taste of what the qualitative behavior of more general solutions may be. But a true qualitative understanding of the nature of general solutions has been achieved only in a neighborhood of the very simplest solutions. The question of the stability of the black hole solutions described above remains unanswered, as do the cosmic censorship conjectures and the nature of the singularities that occur generically in general relativity. Yet these questions are fundamental to the physical interpretation of the theory, and indeed to assessing its very validity.

How likely is it that these questions can ever be answered by rigorous mathematics? Problems concerning the singular behavior of nonlinear hyperbolic partial differential equations are notoriously difficult. The rich geometric structure of the Einstein equations appears at first as a formidable additional complication, but it may also turn out to be a blessing. One can only hope that the Einstein equations will continue to reveal beautiful mathematical structure that answers fundamental questions about our physical world.

# Further Reading

Christodoulou, D. 1999. On the global initial value problem and the issue of singularities. Classical Quantum Gravity 16:A23–A35.   
Hawking, S. W., and G. F. R. Ellis. 1973. The Large Scale Structure of Space-Time. Cambridge Monographs on Mathematical Physics, number 1. Cambridge: Cambridge University Press.   
Penrose, R. 1965. Gravitational collapse and space-time singularities. Physical Review Letters 14:57–59.   
Rendall, A. 2008. Partial Differential Equations in General Relativity. Oxford: Oxford University Press.   
Weyl, H. 1919. Raum, Zeit, Materie. Berlin: Springer. (Also published in English, in 1952, as Space, Time, Matter. New York: Dover.)

# IV.14 Dynamics

# Bodil Branner

# 1 Introduction

Dynamical systems are used to describe the way systems evolve in time, and have their origin in the laws of nature that newton [VI.14] formulated in Principia Mathematica (1687). The associated mathematical discipline, the theory of dynamics, is closely related to many parts of mathematics, in particular analysis, topology, measure theory, and combinatorics. It is also highly influenced and stimulated by problems from the natural sciences, such as celestial mechanics, hydrodynamics, statistical mechanics, meteorology, and other parts of mathematical physics, as well as reaction chemistry, population dynamics, and economics.

Computer simulations and visualizations play an important role in the development of the theory; they have changed our views about what should be considered typical, rather than special and atypical.

There are two main branches of dynamical systems: continuous and discrete. The main focus of this paper will be holomorphic dynamics, which concerns discrete dynamical systems of a special kind. These systems are obtained by taking a holomorphic function [I.3 §5.6] f defined on the complex numbers and applying it repeatedly. An important example is when f is a quadratic polynomial.

# 1.1 Two Basic Examples

It is interesting to note that both types of dynamical system, continuous and discrete, can be well illustrated by examples that date back to Newton.

(i) The N-body problem models the motion in the solar system of the sun and N 1 planets, and does so in terms of differential equations. Each body is represented by a single point, namely its center of mass, and the motion is determined by Newton’s universal law of gravitation—also called the inverse square law. This says that the gravitational force between two bodies is proportional to each of their masses and inversely proportional to the square of the distance between them. Let $r _ { i }$ denote the position vector of the ith body, mi its mass, and g the universal gravitational constant. Then the force on the ith body due to the jth has magnitude $g m _ { i } m _ { j } / \| \pmb { r } _ { j } - \pmb { r } _ { i } \| ^ { 2 }$ , and its direction is along the line from $r _ { i } \mathrm { t o } r _ { j }$ . We can work out the total force on the ith body by adding up all these forces for j  i. Since a unit vector in the direction from ri to rj is $( \pmb { r } _ { j } - \pmb { r } _ { i } ) / | | \pmb { r } _ { j } - \pmb { r } _ { i } |$ , we obtain a force of

$$
g \sum_ {j \neq i} m _ {i} m _ {j} \frac {\boldsymbol {r} _ {j} - \boldsymbol {r} _ {i}}{\| \boldsymbol {r} _ {j} - \boldsymbol {r} _ {i} \| ^ {3}}.
$$

(There is a cube on the bottom rather than a square in order to compensate for the magnitude of $r _ { j } - r _ { i \cdot } )$ A solution to the N-body problem is a set of differentiable vector functions $( r _ { 1 } ( t ) , \ldots , r _ { N } ( t ) )$ ), depending on time t, that satisfy the N differential equations

$$
m _ {i} \boldsymbol {r} _ {i} ^ {\prime \prime} (t) = g \sum_ {j \neq i} m _ {i} m _ {j} \frac {\boldsymbol {r} _ {j} (t) - \boldsymbol {r} _ {i} (t)}{\| \boldsymbol {r} _ {j} (t) - \boldsymbol {r} _ {i} (t) \| ^ {3}},
$$

which result from Newton’s second law, which states that force mass acceleration.

Newton was able to solve the two-body problem explicitly. By neglecting the influence of other planets, he derived the laws formulated by Johannes Kepler, which describe how each planet moves in an elliptic orbit around the sun. However, the jump to $N \ > \ 2$ makes an enormous difference to the complication of the problem: except in very special cases, the system of equations can no longer be solved explicitly (see the three-body problem [V.33]). Nevertheless, Newton’s equations are of great practical importance when it comes to guiding satellites and other space missions.

(ii) newton’s method [II.4 §2.3] for solving equations is quite different and does not involve differential equations. We consider a differentiable function $f$ of one real variable and wish to determine a zero of $f ,$ that is, a solution to the equation $f ( x ) = 0$ . Newton’s idea was to define a new function:

$$
N _ {f} (x) = x - \frac {f (x)}{f ^ {\prime} (x)}.
$$

To put this more geometrically, $N _ { f } ( x )$ is the x-coordinate of the point where the tangent line to the graph $y \ = \ f ( x )$ at the point $( x , f ( x ) )$ crosses the x-axis. (If $f ^ { \prime } ( x ) = 0$ , then this tangent line is horizontal and $N _ { f } ( x )$ is not defined.)

Under many circumstances, if x is close to a zero of $f ,$ then $N _ { f } ( x )$ is significantly closer. Therefore, if we start with some value $x _ { 0 }$ and form the sequence obtained by repeated application of $N _ { f } ,$ that ${ \mathrm { i } } \mathbf { s } ,$ the sequence $x _ { 0 } , x _ { 1 } , x _ { 2 } , \ldots ,$ , where $x _ { 1 } = N _ { f } ( x _ { 0 } ) , x _ { 2 } =$ $N _ { f } ( x _ { 1 } )$ ), and so on, we can expect that this sequence will converge to a zero of $f .$ And this is true: if the initial value $x _ { 0 }$ is sufficiently close to a zero, then the sequence does indeed converge toward that zero, and does so extremely quickly, basically doubling the number of correct digits in each step. This rapid convergence makes Newton’s method very useful for numerical computations.

# 1.2 Continuous Dynamical Systems

We can think of a continuous dynamical system as a system of first-order differential equations, which determine how the system evolves in time. A solution is called an orbit or trajectory, and is parametrized by a number t, which one usually thinks of as time, that takes real values and varies continuously: hence the name “continuous” dynamical system. A periodic orbit of period T is a solution that repeats itself after time $T ,$ , but not earlier.

The differential equation $x ^ { \prime \prime } ( t ) \ = \ - x ( t )$ is of second order, but it is nevertheless a continuous dynamical system because it is equivalent to the system of two first-order differential equations $x _ { 1 } ^ { \prime } ( t ) = x _ { 2 } ( t )$ and $x _ { 2 } ^ { \prime } ( t ) = - x _ { 1 } ( t )$ . In a similar way, the system of differential equations of the N-body problem can be brought into standard form by introducing new variables. The equations are equivalent to a system of 6N first-order differential equations in the variables of the position vectors $\pmb { r } _ { i } = ( x _ { i 1 } , x _ { i 2 } , x _ { i 3 } )$ and the velocity vectors $\pmb { r } _ { i } ^ { \prime } = ( y _ { i 1 } , y _ { i 2 } , y _ { i 3 } )$ . Thus, the N-body problem is a good example of a continuous dynamical system.

In general, if we have a dynamical system consisting of n equations, then we can write the ith equation in the form

$$
x _ {i} ^ {\prime} (t) = f _ {i} (x _ {1} (t), \dots , x _ {n} (t)),
$$

or alternatively we can write all the equations at once in the form ${ \pmb x } ^ { \prime } ( t ) \ = \ { \pmb f } ( { \pmb x } ( t ) )$ , where x(t) is the vector $( x _ { 1 } ( t ) , \ldots , x _ { n } ( t ) )$ and $\pmb { f } = ( f _ { 1 } , \ldots , f _ { n } )$ is a function from $\mathbb { R } ^ { n }$ to $\mathbb { R } ^ { n }$ . Note that f is assumed not to depend on t. If it does, then the system can be brought into standard form by adding the variable $x _ { n + 1 } = t$ and the differential equation $x _ { n + 1 } ^ { \prime } ( t ) = 1$ , which increases the dimension of the system from n to $n + 1$ .

The simplest systems are linear ones, where $f$ is a linear map: that is, $\pmb { f } ( \pmb { x } )$ is given by Ax for some constant n n matrix A. The system above, $x _ { 1 } ^ { \prime } ( t ) =$ $x _ { 2 } ( t )$ and $x _ { 2 } ^ { \prime } ( t ) \ = \ - x _ { 1 } ( t )$ , is an example of a linear system. Most systems, however, including the one for the N-body problem, are nonlinear. If the function $f$ is $" \mathrm { n i c e } ^ { \prime \prime }$ (for instance, differentiable), then uniqueness and existence of solutions are guaranteed for any initial point $\scriptstyle { \pmb x } _ { 0 }$ . That ${ \mathrm { i } } s ,$ there is exactly one solution that passes through the point $\scriptstyle { \pmb x } _ { 0 }$ at time $t = 0 ,$ . For example, in the N-body problem there is exactly one solution for any given set of initial position vectors and initial velocity vectors. It also follows from uniqueness that any pair of orbits must either coincide or be totally disjoint. (Bear in mind that the word “orbit” in this context does not mean the set of positions of a single point mass, but rather the evolution of the vector that represents all the positions and velocities of all the masses.)

Although it is seldom possible to express solutions to nonlinear systems explicitly, we know that they exist, and we call the dynamical system deterministic since solutions are completely determined by their initial conditions. For a given system and given initial conditions it is therefore theoretically possible to predict its entire future evolution.

# 1.3 Discrete Dynamical Systems

A discrete dynamical system is a system that evolves in jumps: “time,” in such a system, is best represented by an integer rather than a real number. A good example is Newton’s method for solving equations. In this instance, the sequence of points we saw earlier, $x _ { 0 } , x _ { 1 } , \ldots , x _ { k } , \ldots $ where $x _ { k } = N _ { f } ( x _ { k - 1 } )$ , is called the orbit of $x _ { 0 }$ . We say that it is obtained by iteration of the function $N _ { f }$ , i.e., by repeated application of the function.

This idea can easily be generalized to other mappings $F : X \to X ,$ , where X could be the real axis, an interval in the real axis, the plane, a subset of the plane, or some more complicated space. The important thing is that the output $F ( x )$ of any input x can be used as the next input. This guarantees that the orbit of any $x _ { 0 }$ in X is defined for all future times. That is, we can define a sequence, $x _ { 0 } , x _ { 1 } , \ldots , x _ { k } , \ldots ,$ where $x _ { k } ~ = ~ F ( x _ { k - 1 } )$ for every k. If the function F has an inverse $F ^ { - 1 }$ , then we can iterate both forwards and backwards and obtain the full orbit of $x _ { 0 }$ as the bi-infinite sequence $\ldots , x _ { - 2 } , x _ { - 1 } , x _ { 0 } , x _ { 1 } , x _ { 2 } , \ldots ,$ , where $x _ { k } = F ( x _ { k - 1 } )$ and, equivalently, $x _ { k - 1 } = F ^ { - 1 } ( x _ { k } )$ , for all integer values.

The orbit of $x _ { 0 }$ is periodic of period k if it repeats itself after time k, but not earlier, i.e., if $x _ { k } = x _ { 0 }$ , but $x _ { j } ~ \neq ~ x _ { 0 }$ for $j = 1 , \ldots , k - 1$ . The orbit is called preperiodic if it is eventually periodic, in other words if there exist $\ell \geqslant 1$ and $k \geqslant 1$ such that $x _ { \ell }$ is periodic of period $k ,$ but none of the $x _ { j }$ for $0 \leqslant j < \ell$ are periodic. The notion of pre-periodicity has no counterpart in continuous dynamics.

A discrete dynamical system is deterministic, since the orbit of any given initial point $x _ { 0 }$ is completely determined once you know $x _ { 0 }$ .

# 1.4 Stability

The modern theory of dynamics was greatly influenced by the work of poincaré [VI.61], and in particular by his prize-winning memoir on the three-body problem, succeeded by three more elaborate volumes on celestial mechanics, all from the late nineteenth century. The memoir was written in response to a competition where one of the proposed problems concerned stability of the solar system. Poincaré introduced the so-called restricted three-body problem, where the third body is assumed to have an infinitely small mass: it does not influence the motion of the other two bodies but it is influenced by them. Poincaré’s work became the prelude to topological dynamics, which focuses on topological properties of solutions to dynamical systems and takes a qualitative approach to them.

Of special interest is the long-term behavior of a system. A periodic orbit is called stable if all orbits through points sufficiently close to it stay close to it at all future times. It is called asymptotically stable if all sufficiently close orbits approach it as time tends to infinity. Let us illustrate this by two linear examples in discrete dynamics. For the real function $F ( x ) = - x ,$ , all points have a periodic orbit: 0 has period 1 and all other x have period 2. Every orbit is stable, but none is asymptotically stable. The real function $\begin{array} { r } { G ( x ) = \frac { 1 } { 2 } x } \end{array}$ has only one periodic orbit, namely 0. Since $G ( 0 ) = 0$ , this orbit has period 1, and we call it a fixed point. If you take any number and repeatedly divide it 2, then the resulting sequence will approach 0, so the fixed point 0 is asymptotically stable.

One of the methods introduced by Poincaré during his study of the three-body problem was a reduction from a continuous dynamical system, in dimension n, say, to an associated discrete dynamical system, a mapping in dimension $n - 1 .$ . The idea is as follows. Suppose we have a periodic orbit of period $T > 0$ in some continuous system. Choose a point $\scriptstyle { \pmb x } _ { 0 }$ on the orbit and a hypersurface Σ through $^ { x _ { 0 } , }$ for instance part of a hyperplane, such that the orbit cuts through Σ at $\scriptstyle { \pmb { x } } _ { 0 } .$ . For any point in Σ that is sufficiently close to $\scriptstyle { \mathbf { x } } _ { 0 } ,$ , one can follow its orbit around and see where it next intersects Σ. This defines a transformation, known as the Poincaré map, which takes the original point to the next point of intersection of its orbit with Σ. It follows from the fact that dynamical systems have unique solutions that every Poincaré map is injective in the neighborhood of $\scriptstyle { \pmb x } _ { 0 }$ (within Σ) for which the Poincaré map is defined. One can perform both forwards and backwards iterations. Note that the periodic orbit of x0 in the continuous system is stable (respectively, asymptotically stable) exactly when the fixed point $\scriptstyle { \pmb x } _ { 0 }$ of the Poincaré map in the discrete system is stable (respectively, asymptotically stable).

# 1.5 Chaotic Behavior

The notion of chaotic dynamics arose in the 1970s. It has been used in different settings, and there is no single definition that covers all uses of the term. However, the property that best characterizes chaos is the phenomenon of sensitive dependence on initial conditions.

Poincaré was the first to observe sensitivity to initial conditions in his treatment of the three-body problem.

Instead of describing his observations let us look at a much simpler example from discrete dynamics. Take as a dynamical space X the half-open unit interval [0, 1), and let F be the function that doubles a number and reduces it modulo 1. That is, $F ( x ) = 2 x$ when $\textstyle 0 \leqslant x < { \frac { 1 } { 2 } }$ and $F ( x ) = 2 x - 1$ when $\begin{array} { r } { \frac { 1 } { 2 } \leqslant x < 1 } \end{array}$ . Let x0 be a number in X and let its iterates be $x _ { 1 } ~ = ~ F ( x _ { 0 } )$ , $x _ { 2 } = F ( x _ { 1 } )$ , and so on. Then $x _ { k }$ is the fractional part of $2 ^ { k } x _ { 0 }$ . (The fractional part of a real number t is what you get when you subtract the largest integer less than t.)

A good way to understand the behavior of the sequence $x _ { 0 } , x _ { 1 } , x _ { 2 } , \ldots$ of iterates is to consider the binary expansion of $x _ { 0 } .$ Suppose, for example, that this begins 0.110100010100111 . . . . To double a number when it is written in binary, all you have to do is shift every digit to the left (just as one does in the decimal system when multiplying by 10). So $2 x _ { 0 }$ will have a binary expansion that begins 1.10100010100111 . . . . To obtain $F ( x _ { 0 } )$ , we have to take the fractional part of this, which we do by subtracting the initial 1. This gives us $x _ { 1 } ~ = ~ 0 . 1 0 1 0 0 0 1 0 1 0 0 1 1 1 \ldots .$ . Repeating the process we find that $x _ { 2 } = 0 . 0 1 0 0 0 1 0 1 0 0 1 1 1 \ldots , x _ { 3 } =$ 0.100010100111 . . . , and so on. (Notice that when we calculated $x _ { 3 }$ from x2 there was no need to subtract 1, since the first digit after the “decimal point” was a 0.) Now consider a different choice of initial number, $x _ { 0 } ^ { \prime } =$ $0 . 1 1 0 1 0 0 0 1 0 1 1 0 1 1 0 \ldots$ . The first nine digits after the decimal point are the same as the first nine digits of x0, so $ { \boldsymbol { { x } } } _ { 0 } ^ { \prime }$ is very close to $x _ { 0 } .$ . However, if we apply F ten times to $x _ { 0 }$ and $\boldsymbol { x } _ { 0 } ^ { \prime } ,$ , then their respective eleventh digits have shifted leftwards and become the first digits of $x _ { 1 0 } = 0 . 0 0 1 1 1 \dots$ . and $x _ { 1 0 } ^ { \prime } = 0 . 1 0 1 1 0 \ldots$ . These two numbers differ by almost ${ \frac { 1 } { 2 } } ,$ , so they are not at all close.

In general, if we know $x _ { 0 }$ to an accuracy of k binary digits and no more, then after k iterations of the map $F$ we have lost all information: $x _ { k }$ could lie anywhere in the interval [0, 1). Therefore, even though the system is deterministic, it is impossible to predict its long-term behavior without knowing $x _ { 0 }$ with perfect accuracy.

This is true in general: it is impossible to make longterm predictions in any part of a dynamical system that shows sensitivity to initial conditions unless the initial conditions are known exactly. In practical applications this is never the case. For instance, when applying a mathematical model to perform weather forecasts, one does not know the initial conditions exactly, and this is why reliable long-term forecasting is impossible.

Sensitivity is also important in the notion of so-called strange attractors. A set A is called an attractor if all orbits that start in A stay in A and if all orbits through nearby points get closer and closer to A. In continuous systems, some simple sets that can be attractors are equilibrium points, periodic orbits (limit cycles), and surfaces such as a torus. In contrast to these examples, strange attractors have both complicated geometry and complicated dynamics: the geometry is fractal and the dynamics sensitive. We shall see examples of fractals later on.

The best-known strange attractor is the Lorenz attractor. In the early 1960s, the meteorologist Edward N. Lorenz studied a three-dimensional continuous dynamical system that gave a simplified model of heat flow. While doing $\mathbf { s o } ,$ he noticed that if he restarted his computer with its initial conditions chosen as the output of an earlier calculation, then the trajectory started to diverge from the one he had previously observed. The explanation he found was that the computer used more precision in its internal calculations than it showed in its output. For this reason, it was not immediately apparent that the initial conditions were in fact very slightly different from before. Because the system was sensitive, this tiny difference eventually made a much bigger difference. He coined the poetic phrase “the butterfly $\mathrm { e f f e c t } ^ { \prime \prime }$ to describe this phenomenon, suggesting that a small disturbance such as a butterfly flickering its wings could in time have a dramatic effect on the longterm evolution of the weather and trigger a tornado thousands of miles away. Computer simulations of the Lorenz system indicate that solutions are attracted to a complicated set that “looks $\mathrm { l i k e } ^ { , \prime }$ a strange attractor. The question of whether it actually was one remained open for a long time. It is not obvious how trustworthy computer simulations are when one is studying sensitive systems, since the computer rounds off the numbers in each step. In 1998 Warwick Tucker gave a computer-assisted proof that the Lorenz attractor is in fact a strange attractor. He used interval arithmetic, where numbers are represented by intervals and estimates can be made precise.

For topological reasons, sensitivity to initial conditions is possible for continuous dynamical systems only when the dimension is at least 3. For discrete systems where the map F is injective, the dimension must be at least 2. However, for noninjective mappings, sensitivity can occur for one-dimensional systems, as we saw with the example given earlier. This is one of the reasons that discrete one-dimensional dynamical systems have been intensively studied.

# 1.6 Structural Stability

Two dynamical systems are said to be topologically equivalent if there is a homeomorphism (a continuous map with continuous inverse) that maps the orbits of one system onto the orbits of the other, and vice versa. Roughly speaking, this means that there is a continuous change of variables that turns one system into the other.

As an example, consider the discrete dynamical system given by the real quadratic polynomial $F ( x ) \ =$ $4 x ( 1 - x )$ . Suppose we were to make the substitution $y \ = \ - 4 x \ + \ 2$ . How could we describe the system in terms of $y \smash { ? }$ Well, if we apply F, then we change x to $4 x ( 1 - x )$ , which means that $y = - 4 x + 2$ changes F(x) $\mathrm { t o } - 4 F ( x ) + 2 = - 1 6 x ( 1 - x ) + 2 $ . But

$$
\begin{array}{l} - 1 6 x (1 - x) + 2 = 1 6 x ^ {2} - 1 6 x + 2 \\ = (- 4 x + 2) ^ {2} - 2 \\ = y ^ {2} - 2. \\ \end{array}
$$

Therefore, the effect of applying the polynomial function F to x is to apply a different polynomial function to y, namely $Q ( y ) = y ^ { 2 }$ 2. Since the change of variables from x to −4x + 2 is continuous and invertible, one says that the functions F and Q are conjugate.

Because F and Q are conjugate, the orbit of any x0 under F becomes, after the change of variables, the orbit of the corresponding point $y _ { 0 } ~ = ~ - 4 x _ { 0 } + 2$ under Q. That is, for every k we have $y _ { k } = - 4 x _ { k } + 2$ . The two systems are topologically equivalent: if you want to understand the dynamics of one of them, you can if you study the other, since its dynamics will be qualitatively the same.

For continuous dynamical systems the notion of equivalence is slightly looser in that we allow a homeomorphism between two topologically equivalent systems to map one orbit onto another without respecting the exact time evolution, but for discrete dynamical systems we must demand that the time evolution is respected as in the example above: in other words, we insist on conjugacy.

The term dynamical system was coined by Stephen Smale in the 1960s and has taken off since then. Smale evolved the theory of robust systems, also named structurally stable systems, a notion that was introduced in the 1930s by Alexander A. Andronov and Lev S. Pontryagin. A dynamical system is called structurally stable if all systems sufficiently close to it, belonging to some specified family of systems, are in fact topologically equivalent to it. We say that they all have the same qualitative behavior. An example of the kind of family one might consider is the set of all real quadratic polynomials of the form $x ^ { 2 } + a$ . This family is parametrized by $^ { a , }$ and the systems close to a given polynomial $x ^ { 2 } + a _ { 0 }$ are all the polynomials $x ^ { 2 } + a$ for which a is close to $\scriptstyle a _ { 0 } .$ We shall return to the question of structural stability when we discuss holomorphic dynamics later.

If a family of dynamical systems parametrized by a variable a is not structurally stable, it may still be that the system with parameter $_ { a _ { 0 } }$ is topologically equivalent to all systems with parameter a in some region that contains $\scriptstyle a _ { 0 } .$ . A major goal of research into dynamics is to understand not just the qualitative structure of each system in the family, but also the structure of the parameter space, that is, how it is divided up into such regions of stability. The boundaries that separate these regions form what is called the bifurcation set: if $_ { a _ { 0 } }$ belongs to this set, then there will be parameters a arbitrarily close to $_ { a _ { 0 } }$ for which the corresponding system has a different qualitative behavior.

A description and classification of structurally stable systems and a classification of possible bifurcations is not within reach for general dynamical systems. However, one of the success stories in the subject, holomorphic dynamics, studies a special class of dynamical systems for which many of these goals have been attained. It is time to turn our attention to this class.

# 2 Holomorphic Dynamics

Holomorphic dynamics is the study of discrete dynamical systems where the map to be iterated is a holomorphic function [I.3 §5.6] of the complex numbers [I.3 §1.5]. Complex numbers are typically denoted by z. In this article, we shall consider iterations of complex polynomials and rational functions (that is, functions like $( z ^ { 2 } + 1 ) / ( z ^ { 3 } + 1 )$ that are ratios of polynomials), but much of what we shall say about them is true for more general holomorphic functions, such as exponential [III.25] and trigonometric [III.92] functions.

Whenever one restricts attention to a special kind of dynamical system, there will be tools that are specially adapted to that situation. In holomorphic dynamics these tools come from complex analysis. When we concentrate on rational functions, there are more special tools, and if we restrict further to polynomials, then there are yet others, as we shall see.

Why might one be interested in iterating rational functions? One answer arose in 1879, when cayley [VI.46] had the idea of trying to find roots of complex polynomials by extending Newton’s method, which we discussed in the introduction, from real numbers to complex numbers. Given any polynomial P, the corresponding Newton function $N _ { P }$ is a rational function, given by the formula

$$
N _ {P} (z) = z - \frac {P (z)}{P ^ {\prime} (z)} = \frac {z P ^ {\prime} (z) - P (z)}{P ^ {\prime} (z)}.
$$

To apply Newton’s method, one iterates this rational function.

The study of the iteration of rational functions flourished at the beginning of the twentieth century, thanks in particular to work of Pierre Fatou and Gaston Julia (who independently obtained many of the same results). Part of their work concerned the study of the local behavior of functions in the neighborhoods of a fixed point. But they were also concerned about global dynamical properties and were inspired by the theory of so-called normal families, then recently established by Paul Montel. However, research on holomorphic dynamics almost came to a stop around 1930, because the fractal sets that lay behind the results were so complicated as to be almost beyond imagination. The research came back to life in around 1980 with the vastly extended calculating powers of computers, and in particular the possibility of making sophisticated graphic visualizations of these fractal sets. Since then, holomorphic dynamics has attracted a lot of attention. New techniques continue to be developed and introduced.

To set the scene, let us start by looking at one of the simplest of polynomials, namely $z ^ { 2 }$ .

# 2.1 The Quadratic Polynomial $z ^ { 2 }$

The dynamics of the simplest quadratic polynomial, $Q _ { 0 } ( z ) ~ = ~ z ^ { 2 }$ , plays a fundamental role in the understanding of the dynamics of any quadratic polynomial. Moreover, the dynamical behavior of $Q _ { 0 }$ can be analyzed and understood completely.

一 ${ \mathrm { ~ f ~ } } z = r \mathrm { e } ^ { \mathrm { i } \theta }$ , then $z ^ { 2 } = r ^ { 2 } \mathrm { e } ^ { 2 \mathrm { i } \theta }$ , so squaring a complex number squares its modulus and doubles its argument. Therefore, the unit circle (the set of complex numbers of modulus 1) is mapped by $Q _ { 0 }$ to itself, while a circle of radius $r \gets 1$ is mapped onto a circle closer to the origin, and a circle of radius $r > 1$ is mapped onto a circle farther away.

Let us look more closely at what happens to the unit circle. A typical point in the circle, $\mathrm { e } ^ { \mathrm { i } \theta }$ , can be parametrized by its argument $\theta ,$ which we can take to lie in the interval [0, 2π). When we square this number, we obtain $\mathrm { e } ^ { 2 \mathrm { i } \theta }$ , which is parametrized by the number 2θ if $2 \theta < 2 \pi ,$ but i $: 2 \theta \geqslant 2 \pi$ , then we subtract 2π so that the argument, $2 \theta - 2 \pi$ , still lies in $[ 0 , 2 \pi )$ . This is strongly reminiscent of the dynamical system we considered in section 1.5. In fact, if we replace the argument θ by its modified argument $\theta / 2 \pi ,$ which amounts to writing e2πiθ instead of $\mathrm { e } ^ { \mathrm { i } \theta }$ , then it becomes exactly the same system. Therefore, the behavior of $z ^ { 2 }$ on the unit circle is chaotic.

As for the rest of the complex plane, the origin is an asymptotically stable fixed point, $Q _ { 0 } ( 0 ) = 0 \ :$ . For any point $z _ { 0 }$ inside the unit circle the iterates $z _ { k }$ converge to 0 as k tends to infinity. For any point $z _ { 0 }$ outside the unit circle the distance zk between the iterates zk and the origin tends to infinity as k tends to infinity. The set of initial points $z _ { 0 }$ with bounded orbit is equal to the closed unit disk, i.e., all points for which $| z _ { 0 } | \leqslant 1$ . Its boundary, the unit circle, divides the complex plane into two domains with qualitatively different dynamical behavior.

Some orbits of $Q _ { 0 }$ are periodic. In order to determine which ones, we first notice that the only possibility outside the unit circle is the fixed point at the origin, since all other points, when you repeatedly square them, either get steadily closer and closer to the origin, or get steadily farther and farther away. So now let us look at the unit circle, and consider the point $\mathrm { e } ^ { 2 \pi \mathrm { i } \theta _ { 0 } }$ , with modified argument $\theta _ { 0 }$ . If this point is periodic with period $k ,$ we must have $2 ^ { k } \theta _ { 0 } = \theta _ { 0 } ( { \mathrm { m o d } } 1 ) { \mathrm { : } }$ that is, $( 2 ^ { k } - 1 ) \theta _ { 0 }$ must be an integer. Because of this, it is convenient to parametrize a point on the unit circle by its modified argument. From now on, when we say “the point $\theta , "$ we shall mean the point $\mathrm { e } ^ { 2 \pi \mathrm { i } \theta }$ , and when we say “argument” we shall mean modified argument.

We have just established that the point θ is periodic with period k only if $( 2 ^ { k } - 1 ) \theta$ is an integer. It follows that there is one point of period 1, namely $\theta _ { 0 } ~ = ~ 0$ . There ait, namely nts of period 2, forming. There are six points for ${ \frac { 1 } { 3 } } \mapsto { \frac { 2 } { 3 } } \mapsto { \frac { 1 } { 3 } }$ period 3, forming two orbits, namely $\textstyle { \frac { 1 } { 7 } } \mapsto { \frac { 2 } { 7 } } \mapsto { \frac { 4 } { 7 } } \mapsto { \frac { 1 } { 7 } }$ and ${ \frac { 3 } { 7 } } \mapsto { \frac { 6 } { 7 } } \mapsto { \frac { 5 } { 7 } } \mapsto { \frac { 3 } { 7 } }$ . (At each stage, we double the number we have, and subtract 1 if that is needed to get us back into the interval [0, 1).) The points of period 4 are fractions with denominator 15, but the converse is not true: the fractions $\begin{array} { r } { \frac { 3 } { 1 5 } = \frac { 1 } { 3 } } \end{array}$ and $\textstyle { \frac { 6 } { 1 5 } } = { \frac { 2 } { 3 } }$ have the lower period 2. The periodic points on the unit circle are dense in the unit circle, meaning that arbitrarily close to any point is a periodic point. This follows from the observation that all repeating binary expansions, such as 0.1100011000110001100011000 . . . are periodic, and any finite sequence of 0s and 1s is the start of a repeating sequence. One can, in fact, show that the periodic points on the unit circle are exactly the points whose argument is a fraction $p / q$ in $[ 0 , 1 )$ with q odd. Any fraction with even denominator can be written in the form $p / ( 2 ^ { \ell } q )$ for some odd number q. After  iterations, such a fraction will land on a periodic point, so the initial point is pre-periodic. Points with rational argument in [0, 1) have a finite orbit, while points with irrational argument have an infinite orbit. The reason for taking modified arguments is now justified: the behavior of the dynamics depends on whether $\theta _ { 0 }$ is rational or irrational.

When $\theta _ { 0 }$ is irrational its orbit may or may not be dense in $[ 0 , 1 )$ . This is another fact that is easy to see if one considers binary expansions. For instance, a very special example of a $\theta _ { 0 }$ with a dense orbit is given by the binary expansion

$$
\theta_ {0} = 0. 0 1 0 0 0 1 1 0 1 1 0 0 0 0 0 1 0 1 0 0 1 1 1 0 0 1 0 1 1 1 0 1 1 1 \dots ,
$$

where one obtains this expansion by simply listing all finite binary sequences in turn: first the blocks of length one, 0 and 1, then the blocks of length two, 00, 01, 10, and 11, and so on. When we iterate, this binary expansion shifts to the left and all possible finite sequences appear at some time or another at the beginning of some iterate $\theta _ { k }$ .

# 2.2 Characterization of Periodic Points

Let $z _ { 0 }$ be a fixed point of a holomorphic map F. How do the iterates of points near $z _ { 0 }$ behave? The answer depends crucially on a number $\rho ,$ , called the multiplier of the fixed point, which is defined to be $F ^ { \prime } ( z _ { 0 } )$ . To see why this is relevant, notice that if z is very close to $z _ { \mathrm { 0 } } ,$ then $F ( z )$ is, to a first-order approximation, equal to $F ( z _ { 0 } ) + F ^ { \prime } ( z _ { 0 } ) ( z - z _ { 0 } ) = z _ { 0 } + \rho ( z - z _ { 0 } )$ . Thus, when you apply F to a point near $z _ { 0 } ,$ its difference from $z _ { 0 }$ approximately multiplies by $\rho . \mathrm { ~ I f ~ } | \rho | < 1$ , then nearby points will get closer to $z _ { 0 } ,$ , in which case $z _ { 0 }$ is called an attracting fixed point. If $\rho = 0$ , then this happens very quickly and $z _ { 0 }$ is called super-attracting. If $| \rho | > 1$ , then nearby points get farther away and $z _ { 0 }$ is called repelling. Finally, if $| \rho | = 1$ , then one says that z0 is indifferent.

If $z _ { 0 }$ is indifferent, then its multiplier will take the form $\rho = \mathrm { e } ^ { 2 \pi \mathrm { i } \theta }$ , and near $z _ { 0 }$ the map F will be approximately a rotation about $z _ { 0 }$ by an angle of 2πθ. The behavior of the system depends very much on the precise value of θ. We call the fixed point rationally or irrationally indifferent if θ is rational or irrational, respectively. The dynamics is not yet completely understood in all irrational cases.

A periodic point $z _ { 0 }$ of period k will be a fixed point of the kth iterate $F ^ { k } = F \circ \cdots \circ F$ of F. For this reason we define its multiplier by $\rho = ( F ^ { k } ) ^ { \prime } ( z _ { 0 } ) $ ). It follows from the chain rule that

$$
\left(F ^ {k}\right) ^ {\prime} \left(z _ {0}\right) = \prod_ {j = 0} ^ {k - 1} F ^ {\prime} \left(z _ {j}\right)
$$

and therefore that the derivative of $F ^ { k }$ is the same at all points of the periodic orbit. This formula also implies that a super-attracting periodic orbit must contain a critical point (that is, a point where the derivative of F is zero): if $( F ^ { k } ) ^ { \prime } ( z _ { 0 } ) = 0$ , then at least one $F ^ { \prime } ( z _ { j } )$ must be 0.

Note that 0 is a super-attracting fixed point of $Q _ { 0 }$ , and that any periodic orbit of $Q _ { 0 }$ of period k on the unit circle has multiplier $2 ^ { k }$ . All periodic orbits on the unit circle are therefore repelling.

# 2.3 A One-Parameter Family of Quadratic Polynomials

The quadratic polynomial $Q _ { 0 }$ sits at the center of the one-parameter family of quadratic polynomials of the form $Q _ { c } ( z ) = z ^ { 2 } + c$ . (We considered this family earlier, but then z and c were real rather than complex.) For each fixed complex number c we are interested in the dynamics of the polynomial $Q _ { c }$ c under iteration. The reason we do not need to study more general quadratic polynomials is that they can be brought into this form by a simple substitution $w \ = \ a z \ + \ b$ , similar to the substitution in the real example in section 1.6. For any given quadratic polynomial P we can find exactly one substitution $w = a z + b$ and one c such that

$$
a (P (z)) + b = (a z + b) ^ {2} + c \quad \text { for   all } z.
$$

Therefore, if we understand the dynamics of the polynomials $\boldsymbol { Q _ { c } } ,$ , then we understand the dynamics of all quadratic polynomials.

There are other representative families of quadratic polynomials that can be useful. One example is the family $F _ { \lambda } ( z ) = \lambda z + z ^ { 2 }$ . The substitution $w = z + { \frac { 1 } { 2 } } \lambda$ changes $F _ { \lambda }$ into $Q _ { c }$ , where $\begin{array} { r } { c = \frac { 1 } { 2 } \lambda - \frac { 1 } { 4 } \lambda ^ { 2 } } \end{array}$ . We shall return to the expression of c in terms of λ later on. In the family of polynomials $\displaystyle Q _ { c } ,$ the parameter $c = Q _ { c } ( 0 )$ coincides with the only critical value of $Q _ { c }$ in the plane:

![](images/af6cccc23f7f5df0a6e4f598f2bbe4428c5de58446c507081664c06432948f4d.jpg)

<details>
<summary>text_image</summary>

N
Z
R
z
iR
S
C
</details>

Figure 1 The Riemann sphere.

as we shall see later, critical orbits play an essential role in the analysis of the global dynamics. In the family of polynomials $F _ { \lambda }$ the parameter λ is equal to the multiplier of the fixed point at the origin of $F _ { \lambda } ,$ which sometimes makes this family more convenient.

# 2.4 The Riemann Sphere

To understand further the dynamics of polynomials it is best to regard them as a special case of rational functions. Since a rational function can sometimes be infinite, the natural space to consider is not the complex plane C but the extended complex plane, which is the complex plane together with the point “∞.” This space is denoted ${ \hat { \mathbb { C } } } = \mathbb { C } \cup \{ \infty \}$ . A geometrical picture (see figure 1) is obtained by identifying the extended complex plane with the Riemann sphere. This is simply the unit sphere $\{ ( x _ { 1 } , x _ { 2 } , x _ { 3 } ) : x _ { 1 } ^ { 2 } + x _ { 2 } ^ { 2 } + x _ { 3 } ^ { 2 } = 1 \}$ in three-dimensional space. Given a number z in the complex plane, the straight line joining z to the north pole $N = ( 0 , 0 , 1 )$ ) intersects this sphere in exactly one place (apart from N itself). This place is the point in the sphere that is associated with z. Notice that the bigger z is, the closer the associated point is to N. We therefore regard N as corresponding to the point ∞.

Let us now think of $Q _ { 0 } ( z ) = z ^ { 2 }$ as a function from Cˆ to Cˆ. We have seen that 0 is a super-attracting fixed point of $Q _ { 0 }$ . What about $\infty ,$ , which is a fixed point as well? The classification we gave in terms of multipliers does not work at , but a standard trick in this situation is to “move” to 0. If one wishes to understand the behavior of a function $f$ with a fixed point at , one can look instead at the function $g ( z ) ~ = ~ 1 / f ( 1 / z )$ , which has a fixed point at 0 (since $1 / f ( 1 / 0 ) ~ = ~ 1 / f ( \infty ) ~ =$ $1 / \infty \ = \ 0 )$ . When $f ( z ) ~ = ~ z ^ { 2 } , ~ g ( z )$ is also $z ^ { 2 } ,$ , so  is also a super-attracting fixed point of $Q _ { 0 }$ .

![](images/c4955e35f2ece18fa29ee1ad140d9072eab1e235a0613f796251072b9813b0cf.jpg)

<details>
<summary>natural_image</summary>

Abstract fractal pattern with black and gray polygonal shapes (no text or symbols)
</details>

Figure 2 The Douady rabbit. The filled Julia set of $Q _ { c _ { 0 } }$ where c0 is the one root of the polynomial $( c ^ { 2 } + c ) ^ { 2 } + { \bar { c } }$ that has positive imaginary part. This corresponds to one of the three possible c values for which the critical orbit $0 \mapsto c \mapsto c ^ { 2 } + c \mapsto ( c ^ { 2 } + c ) ^ { 2 } + c = 0$ is periodic of period 3. The critical orbit is marked with three white dots inside the filled Julia set: 0 in the black, $c _ { 0 }$ in the light gray, and $c _ { 0 } ^ { 2 } + c _ { 0 }$ in the gray. The corresponding three attracting basins of $Q _ { c _ { 0 } } ^ { 3 }$ Q3c are marked in black, light gray, and gray, respectively. The Julia set is the common boundary of the black, light gray, and gray basins of attraction as well as of $A _ { c _ { 0 } } ( \infty )$ ).

In general, if P is any nonconstant polynomial, then it is natural to define $P ( \infty )$ to be . Applying the above trick, we obtain a rational function. For example, if $P ( z ) = z ^ { 2 } + 1$ , then $1 / P ( 1 / z ) = z ^ { 2 } / ( z ^ { 2 } + 1 )$ . If P has degree at least 2, then is a super-attracting fixed point.

The connection between $\hat { \mathbb { C } }$ and rational functions is expressed by the following fact: a function $F : \hat { \mathbb { C } } \to \hat { \mathbb { C } }$ is holomorphic everywhere (with suitable definitions at ) if and only if it is a rational function. This is not obvious, but is typically proved in a first course in complex analysis. Among the rational functions, the polynomials are the ones for which $F ( \infty ) = \infty =$ $F ^ { - 1 } ( \infty )$ .

A polynomial P of degree d has $d - 1$ critical points in the plane (not including ). These are the roots of the derivative $P ^ { \prime }$ , counted with multiplicity. The critical point at ∞ has multiplicity $d - 1$ , as can again be seen by looking at the map $1 / P ( 1 / z )$ . In particular, quadratic polynomials have exactly one critical point in the plane. The degree of a rational function $P / Q$ (where P and Q have no common roots) is defined to be the maximal degree of the polynomials P and Q. A rational function of degree d has $2 d - 2$ critical points in $\hat { \mathbb { C } } ,$ as we have just seen for polynomials.

# 2.5 Julia Sets of Polynomials

It can be shown that the only invertible holomorphic maps from $\mathbb { C } \operatorname { t o } \mathbb { C }$ are polynomials of degree 1, that is, functions of the form $a z + b$ with $\alpha \neq 0 .$ . The dynamical behavior of these maps is easy to analyze, simple, and hence not interesting.

From now on, therefore, we shall consider only polynomials P of degree at least 2. For all such polynomials, ∞ is a super-attracting fixed point, from which it follows that the plane is split into two disjoint sets with qualitatively different dynamics, one consisting of points that are attracted to and the other consisting of points that are not. The attracting basin of $\infty ,$ denoted by $A _ { P } ( \infty )$ , consists of all initial points z such that $P ^ { k } ( z ) \to \infty$ as $k  \infty$ . (Here, $P ^ { k } ( z )$ stands for the result of applying P to z k times.) The complement of $A _ { P } ( \infty )$ is called the filled Julia set, and is denoted by $K _ { P }$ . It can be defined as the set of all points z such that the sequence $z , P ( z ) , P ^ { 2 } ( z ) , P ^ { 3 } ( z ) , . .$ . is bounded. (It is not hard to show that sequences of this kind either tend to or are bounded.)

The attracting basin of  is an open set and the filled Julia set is a closed, bounded set (i.e., a compact set [III.9]). The attracting basin of ∞ is always connected. For this reason the boundary of $K _ { P }$ is equal to the boundary of $A _ { P } ( \infty )$ . The common boundary is called the Julia set of P and is denoted by $J _ { P } ,$ . The three sets $K _ { P } , A _ { P } ( \infty )$ , and $J _ { P }$ are completely invariant, i.e., $P ( K _ { P } ) = K _ { P } = P ^ { - 1 } ( K _ { P } )$ , and so on. If we replace P by any iterate $P ^ { k } ,$ , then the filled Julia set, the attracting basin of $\infty ,$ , and the Julia set of $P ^ { k }$ are the same sets as those of P.

For the polynomial $Q _ { 0 } ,$ , we showed earlier that the filled Julia set is the closed unit disk, $\{ z : | z | \leqslant 1 \} ;$ ; the attracting basin of  is its complement, $\{ z : | z | > 1 \}$ ; and the Julia set is the unit circle, $\{ z : | z | = 1 \}$ .

The name “filled Julia set” refers to the fact that $K _ { P }$ is equal to $J _ { P }$ with all its holes (or, more formally, the bounded components of its complement) filled in. The complement of the Julia set is called the Fatou set and any connected component of it is called a Fatou component.

Figures 2–6 show different examples of Julia sets of quadratic polynomials $Q _ { c }$ . For simplicity we set $K _ { Q _ { c } } = K _ { c } , A _ { Q _ { c } } ( \infty ) = A _ { c } ( \infty )$ , and $J _ { Q _ { c } } = J _ { c }$ . Note that all Julia sets $J _ { c }$ are symmetric around 0, owing to the symmetry in the formula: $Q _ { c } ( - z ) = Q _ { c } ( z )$ , which implies that if a point z belongs to $J _ { c }$ , then so does z.

![](images/5c839603a651b4f545dcf55597d63feed0170b55017610280d2e7f7dd488b58b.jpg)

<details>
<summary>natural_image</summary>

Decorative cloud-like shape with no text or symbols
</details>

Figure 3 The Julia set of $Q _ { 1 / 4 }$ . Every point inside the Julia set (including the critical point 0) is attracted (under repeated applications of $Q _ { 1 / 4 } )$ to the rationally indifferent fixed point 1 with multiplier $\rho = 1$ , which belongs to $J _ { 1 / 4 }$ .

![](images/834f774141201fe79dead3812e6be58934d4c4c1649c05259e4c8ef8eca93ab0.jpg)

<details>
<summary>natural_image</summary>

Fractal pattern with concentric circles and irregular star-like shapes (no text or symbols)
</details>

Figure 4 The Julia set of $Q _ { c }$ with a so-called Siegel disk around an irrationally indifferent fixed point of multiplier $\rho ~ = ~ \mathrm { e } ^ { 2 \pi \mathrm { i } ( \sqrt { 5 } - 1 ) / 2 }$ . The corresponding c-value is equal to ${ \textstyle { \frac { 1 } { 2 } } } \rho - { \textstyle { \frac { 1 } { 4 } } } \rho ^ { 2 } .$ In the Siegel disk, the Fatou component containing the fixed point, the action of $Q _ { c }$ can, after a suitable change of variables, be expressed as $w \mapsto \rho w$ . The fixed point is marked and so are some orbits of points in its vicinity. The critical orbit is dense in the boundary of the Siegel disk.

# 2.6 Properties of Julia Sets

In this section we shall list several common properties of Julia sets. The proofs of these, which are beyond the scope of this article, mostly depend on the theory of normal families.

The Julia set is the set of points for which the system displays sensitivity to initial conditions, i.e., the chaotic subset of the dynamical system.

The repelling orbits belong to the Julia set and form a dense subset of the set. That is, any point in the Julia set can be approximated arbitrarily well by a repelling point. This is the definition originally used by Julia. (Of course, the name “Julia set” was used only later.)   
• For any point z in the Julia set, the set of iterated preimages $\textstyle \bigcup _ { k = 1 } ^ { \infty } F ^ { - k } ( z )$ forms a dense subset of the Julia set. This property is used when one is making computer pictures of Julia sets.   
In fact, for any point z in Cˆ (with at most one or two exceptions), the closure of the set of iterated preimages contains the Julia set.   
• For any point z in the Julia set and any neighborhood $U _ { z }$ of $z ,$ the iterated images $F ^ { k } ( U _ { z } )$ cover all of Cˆ except at most one or two exceptional points. This property demonstrates an extreme sensitivity to initial conditions.   
If Ω is a union of Fatou components that is completely invariant (that is, ${ \cal F } ( \Omega ) ~ = ~ \Omega ~ = ~ { \cal F } ^ { - 1 } ( \Omega ) )$ , then the boundary of Ω coincides with the Julia set. This justifies the definition of the Julia set of a polynomial as the boundary of the attracting basin of . Compare also with figure ${ } ^ { 2 , }$ where the attracting basins of $Q _ { c _ { 0 } } ^ { 3 }$ and $A _ { c _ { 0 } } ( \infty )$ ) are examples of such completely invariant sets.   
The Julia set is either connected or consists of uncountably many connected components. An example of the latter is shown in figure 6.   
The Julia set is typically a fractal: when one zooms in on it, one finds that the complication of the set is repeated at all scales. It is also self-similar, in the following sense: for any noncritical point z in the Julia set, any sufficiently small neighborhood $U _ { z }$ of z is mapped bijectively onto $F ( U _ { z } )$ , a neighborhood of $F ( z )$ . The Julia set in $U _ { z }$ and the Julia set in $F ( U _ { z } )$ look alike.

All but the last two properties can easily be verified in the example $Q _ { 0 }$ . In this case the exceptional points are 0 and .

# 2.7 Böttcher Maps and Potentials

# 2.7.1 Böttcher Maps

Consider the quadratic polynomial $Q _ { - 2 } ( z ) = z ^ { 2 } - 2$ . If z belongs to the interval $[ - 2 , 2 ]$ , then $z ^ { 2 }$ belongs to the interval [0, 4], so $Q _ { - 2 } ( z )$ also belongs to the interval [ 2, 2]. It follows that this interval is contained in the filled Julia set $K _ { - 2 }$ .

![](images/6e67cffdfe50d14f614ad6bfca0c0f3cfe9c1ed990f4d499f3ce28466bf2b300.jpg)

<details>
<summary>text_image</summary>

(a)
φ₋₂ ψ₋₂
</details>

![](images/862306bd5559419ecae5e53d977e2cc538155c126fc4bb7e12ed07b34f9cb838.jpg)

<details>
<summary>natural_image</summary>

Pure geometric diagram with concentric ellipses and radial lines, no text or symbols present
</details>

Figure 5 (a) Some equipotentials and external rays $\mathcal { R } _ { 0 } ( \theta )$ of $Q _ { 0 }$ in $A _ { 0 } ( \infty )$ , the set of complex numbers of modulus greater than 1. (b) The corresponding equipotentials and external rays $\mathcal { R } _ { - 2 } ( \theta )$ of $Q _ { - 2 }$ in $A _ { - 2 } ( \infty )$ , the set of complex numbers not in $K _ { - 2 } = J _ { - 2 } = [ - 2 , 2 ]$ . The external rays that are drawn have arguments $\begin{array} { r } { \theta = \frac { 1 } { 1 2 } p , } \end{array}$ , where $p = 0 , 1 , \ldots , 1 1$ .

The polynomial $Q _ { - 2 } ( z )$ is not topologically equivalent to $Q _ { 0 } ( w ) ~ = ~ w ^ { 2 }$ , but when z is big enough, it behaves in a similar way, since 2 is small compared with $z ^ { 2 } .$ . We can express this similarity with an appropriate holomorphic change of variables. Indeed, suppose that $z = w + 1 / w$ . Then when w changes to $w ^ { 2 }$ , z changes to $\textstyle w ^ { 2 } + 1 / w ^ { 2 }$ . But this equals

$$
(w + 1 / w) ^ {2} - 2 = z ^ {2} - 2 = Q _ {- 2} (z).
$$

The reason this does not show that $Q _ { 0 }$ and $Q _ { - 2 }$ are equivalent is that the change of variables cannot be inverted. However, in a suitable region it can. If $z =$ $\boldsymbol { \upsilon } + 1 / \boldsymbol { \upsilon }$ , then $w ^ { 2 } - w z + 1 = 0$ . Solving this quadratic equation we find that $\begin{array} { r } { w = { \frac { 1 } { \gamma } } ( z \pm { \sqrt { z ^ { 2 } - 4 } } ) } \end{array}$ , which leaves us with the problem of which square root to take. It can be shown that for one choice $| \boldsymbol { w } | < 1$ and for the other choice $| \boldsymbol { w } | > 1$ , as long as z does not lie in the interval [−2, 2]. If we always choose the square root for which $| \boldsymbol { w } | > 1$ , then it turns out that the resulting function of z is a continuous function (in fact, holomorphic) from the set C [ 2, 2] of complex numbers not in [ 2, 2] to the set $\{ w : | w | > 1 \}$ of complex numbers of modulus greater than 1.

Once this is established, it follows that the behavior of $Q _ { - 2 }$ on the set $\mathbb { C } \setminus [ - 2 , 2 ]$ is topologically the same as the behavior of $Q _ { 0 }$ on the set $\{ w \ : \ | w | \ > \ 1 \}$ . In particular, points outside $\mathbb { C } \setminus [ - 2 , 2 ]$ have orbits that tend to infinity under iteration by $Q _ { - 2 }$ . Therefore, the attracting basin $A _ { - 2 } ( \infty )$ of $Q _ { - 2 } { \mathrm { ~ i s ~ } } \mathbb { C } \setminus \left[ - 2 , 2 \right]$ , and the filled Julia set $K _ { - 2 }$ and the Julia set $J _ { - 2 }$ are both equal to [ 2, 2].

Let us write $\psi _ { - 2 } ( w )$ for $\boldsymbol { \upsilon } + 1 / \boldsymbol { \upsilon }$ . The function $\psi _ { - 2 } ,$ which we used to change variables, maps circles of radius greater than 1 onto ellipses, and takes radial lines $\mathcal { R } _ { 0 } ( \theta )$ that consists of all complex numbers of some given argument θ and modulus greater than 1 to half-branches of hyperbolas. Since the ratio of $\psi _ { - 2 } ( w )$ to w tends to 1 as $w \to \infty$ , each radial line will be the asymptote of the corresponding hyperbola half-branch (see figure 5).

It turns out that what we have just done for the polynomial $Q _ { - 2 }$ can be done for any quadratic polynomial $Q _ { c }$ . That is, for sufficiently large complex numbers there is a holomorphic function, denoted $\varphi _ { c }$ , called the Böttcher map, that changes variables in such a way that $Q _ { c }$ turns into $Q _ { 0 } ,$ , in the sense that $\varphi _ { c } ( Q _ { c } ( z ) ) =$ $\varphi _ { c } ( z ) ^ { 2 }$ . (The map $\psi _ { - 2 }$ described above is the inverse of the Böttcher map in the case $c \ = \ - 2$ , rather than the map itself.) After the change of variables, the new coordinates are called Böttcher coordinates.

More generally, for all monic polynomials P (i.e., polynomials with leading coefficient 1) there is a unique holomorphic change of variables $\varphi _ { P }$ that converts P into the function $z \mapsto z ^ { d }$ for large enough $z ,$ in the sense that $\varphi _ { P } ( P ( z ) ) = \varphi _ { P } ( z ) ^ { d }$ , and has the property that $( \varphi _ { P } ( z ) / z ) \to 1$ as $z  \infty$ . The inverse of $\varphi _ { P }$ is written $\psi _ { P }$ .

# 2.7.2 Potentials

As we have noted already, if one repeatedly squares a complex number z of modulus greater than 1, then it will escape to infinity. The larger the modulus of $z ,$ the faster the iterates will tend to infinity. If instead of squaring, one applies a monic polynomial P of degree $^ { d , }$ then for large enough z it is again true that the iterates $z , P ( z ) , P ^ { 2 } ( z ) , \ldots$ tend to infinity. It follows from the formula $\begin{array} { l } { \displaystyle \phi _ { P } ( \boldsymbol { P } ( z ) ) ~ = ~ \phi _ { P } ( z ) ^ { d } } \end{array}$ that $\varphi _ { P } ( P ^ { k } ( z ) ) = \varphi _ { P } ( z ) ^ { d ^ { k } }$ . Therefore, the speed at which the iterates tend to infinity depends not on z but on $| \varphi _ { P } ( z ) |$ : the larger the value of $| \varphi _ { P } ( z ) |$ , the faster the convergence. For this reason, the level sets of $| \varphi _ { P } | ,$ that is, sets of the form $\{ z \in \mathbb { C } : | \varphi _ { P } ( z ) | = r \}$ , are important.

For many purposes it is useful to look not at the function $\varphi _ { P }$ itself but at the function $g _ { P } ( z ) = \log | \varphi _ { P } ( z ) |$ . This function is called the potential, or Green’s function. It has the same level sets as $| \varphi _ { P } ( z ) |$ , but has the advantage that it is a harmonic function [IV.24 §5.1].

Clearly, $_ { g _ { P } }$ is defined whenever $\varphi _ { P }$ is defined. But we can in fact extend the definition of $_ { g p }$ to the whole of the attracting basin $A _ { P } ( \infty )$ . Given any z for which the iterates $P ^ { k } ( z )$ tend to infinity, one chooses some k such that $\varphi _ { P } ( P ^ { k } ( z ) )$ is defined and one sets $g _ { P } ( z )$ to be $d ^ { - k }$ log $| \varphi _ { P } ( P ^ { k } ( z ) ) |$ . Notice that $\varphi _ { P } ( P ^ { k + 1 } ( z ) ) =$ $\varphi _ { P } ( P ^ { k } ( z ) ) ^ { d }$ , so log $| \varphi _ { P } ( P ^ { k + 1 } ( z ) ) | = d \log | \varphi _ { P } ( P ^ { k } ( z ) ) | ,$ , from which it is easy to deduce that the value of $d ^ { - k } \log | \varphi _ { P } ( P ^ { k } ( z ) ) |$ does not depend on the choice of k.

The level sets of $_ { g p }$ are called equipotentials. Notice that the equipotential of potential $g _ { P } ( z )$ is mapped by P onto the equipotential of potential $g _ { P } ( P ( z ) ) =$ $d g _ { P } ( z )$ . As we shall see, useful information about the dynamics of the polynomial $P$ can be deduced from information about its equipotentials.

If $\psi _ { P }$ is defined everywhere on the circle $C _ { r }$ of radius $r ,$ for some $r \ > \ 1$ , then it maps it to $\left\{ z \right. :$ $| \varphi _ { P } ( z ) | ~ = ~ r \}$ , which is the equipotential of potential log r . For large enough r , this equipotential is a simple closed curve encircling $K _ { P }$ , and it shrinks as r decreases. It is possible for two parts of this curve to come together so that it forms a figure-of-eight shape and then splits into two, like an amoeba dividing, but this can happen only if the curve crosses a critical point of P . Therefore, if all the critical points of P belong to the filled Julia set $K _ { P }$ (as in the example $Q _ { - 2 }$ where $0 \in K _ { - 2 } = [ - 2 , 2 ] )$ , then it cannot happen. In this case, the Böttcher map $\varphi _ { P }$ can be defined on the whole of the attracting basin $A _ { P } ( \infty )$ , and it is a bijection from $A _ { P } ( \infty )$ to the attracting basin $A _ { 0 } ( \infty ) = \{ w \in \mathbb { C } : | w | > 1 \}$ of the polynomial $z ^ { d }$ . There are equipotentials of potential t for every $t > 0$ and they are all simple closed curves. (Compare with figure 5.) As t approaches 0, the equipotential of potential t, together with its interior, forms a shape that gets closer and closer to the filled Julia set $K _ { P }$ . It follows that $K _ { P }$ is a connected set, as is the Julia set $J _ { P }$ .

![](images/83874a50dcf5f18b597dc52ab9a80430c92b85482639e16514fe39cde3954311.jpg)

<details>
<summary>natural_image</summary>

Abstract symmetrical pattern with floral motifs and a central diagonal line (no text or symbols)
</details>

Figure 6 The Julia set of a quadratic polynomial $Q _ { c }$ for which the critical point 0 escapes to infinity under iteration. The Julia set is totally disconnected. The figure-of-eight-shaped curve with 0 at its intersection point is the equipotential through 0. The simple closed curve surrounding it is the equipotential through the critical value c.

On the other hand, if at least one of the critical points in the plane belongs to $A _ { P } ( \infty )$ , then at a certain point the image of $C _ { r }$ splits into two or more pieces. In particular, the equipotential containing the fastest escaping critical point $( \mathrm { i . e . }$ , the critical point with the highest value of the potential $_ { g _ { P } ) }$ has at least two loops, as is illustrated in figure 6. The inside of each loop is mapped by P onto the inside of the equipotential of the corresponding critical value, which is a simple closed curve (since the potential of the critical value is greater than the potential of any critical point). Inside each loop there must be points from the filled Julia set $K _ { P } ,$ , so this set must be disconnected. The Böttcher map can always be defined on the outside of the equipotential of the fastest escaping critical point and can therefore always be applied to the fastest escaping critical value.

If $Q _ { c }$ is a quadratic polynomial for which 0 escapes to infinity under iteration, then the filled Julia set turns out to be totally disconnected, which means that the connected components of $K _ { c }$ are points. None of these points is isolated: they can all be obtained as limits of sequences of other points of $K _ { c } . \mathrm { A }$ set which is compact, totally disconnected, and with no isolated points is called a cantor set [III.17], since such a set is homeomorphic to Cantor’s middle-thirds set. Note that in this case $K _ { c } = J _ { c }$ . For $Q _ { \ l }$ c we have the following dichotomy: the Julia set $J _ { c }$ is connected if 0 has a bounded orbit, and it is totally disconnected if 0 escapes to infinity under iteration. We shall return to this dichotomy when we come to define the Mandelbrot set later in this article.

# 2.7.3 External Rays of Polynomials with Connected Julia Set

We have just obtained information by looking at the images under $\psi _ { P }$ of circles of radius greater than 1. We can obtain complementary information from the images of radial lines, which cut all these circles at right angles. If the Julia set is connected, then, as we saw in the discussion of potentials, the Böttcher map $\varphi _ { P }$ is a bijection from the attracting basin $A _ { P } ( \infty )$ to the attracting basin of $z ^ { d } ,$ , which is the complement $\{ w : | w | > 1 \}$ of the closed unit disk. As before, let $\mathcal { R } _ { 0 } ( \theta )$ denote the half-line that consists of all complex numbers of argument θ and modulus greater than 1. Because $( \varphi _ { P } ( z ) / z ) \to 1 \mathrm { ~ a s ~ } z \to \infty$ , the image of ${ \mathcal { R } } _ { 0 } ( \theta )$ under $\psi _ { P }$ is a half-infinite curve consisting of points with arguments getting closer and closer to θ. This curve is denoted by $\mathcal { R } _ { P } ( \theta )$ , and is known as the external ray of argument θ of $P .$ Note that $\mathcal { R } _ { 0 } ( \theta )$ is the external ray of argument θ of $z ^ { d }$ .

One can think of equipotentials as contour lines of the potential function, and of external rays as the lines of steepest ascent. Between the two of them, they provide a parametrization of the attracting basin, just as modulus and argument provide a parametrization of $\{ z : | z | > 1 \} :$ if you know the potential at a certain complex number z, and you also know which external ray it lies on, then you know what z is. Moreover, a ray of argument θ is mapped by P onto the ray of argument $d \theta ,$ , just as, when a number z lies on the half-line $\mathcal { R } _ { 0 } ( \theta )$ , then $z ^ { d }$ lies on the half-line $\mathcal { R } _ { 0 } ( d \theta )$ .

We say that an external ray lands if $\psi _ { P } ( r \mathrm { e } ^ { 2 \pi \mathrm { i } \theta } )$ converges to a limit as $r \setminus 1$ . If this happens, then the limit is called the landing point. However, it may happen that the end of the ray oscillates so much that there is a continuum of different limit points. In this case the ray is nonlanding. It can be shown that all rational rays land. Since a rational ray is either periodic or pre-periodic under iteration by P, the landing point of a rational ray must be either a periodic or a pre-periodic point in the Julia set. Much of the structure of the Julia set can be picked up from knowledge about common landing points. In the example illustrated in figure 2, the closures of the three Fatou components containing the critical orbit have one point in common. This point is a repelling fixed point and the common landing point of the rays of argument $\textstyle { \frac { 1 } { 7 } } , \ { \frac { 2 } { 7 } } , \ { \frac { 4 } { 7 } }$ . The rays of argument $\frac { 1 } { 7 }$ and $\frac { 2 } { 7 }$ are adjacent to the Fatou component containing the critical value $c _ { 0 } .$ These two arguments will show up again in the parameter plane and tell us where $c _ { 0 }$ is situated.

# 2.7.4 Local Connectedness

In the example illustrated in figure 5 the inverse of the Böttcher map (the function $\psi _ { - 2 } )$ is defined on the set $\{ w \ : \ | w | \ > \ 1 \}$ of all complex numbers w of modulus greater than 1. However, it can be continuously extended to a function defined on the larger set {w : $| \boldsymbol { w } | \geqslant 1 \}$ . If we use the formula $\psi _ { - 2 } ( w ) = w + 1 / w$ , then we have $\psi _ { - 2 } ( { \mathrm e } ^ { 2 \pi \mathrm { i } \theta } ) = 2 \cos ( 2 \pi \theta ) $ , which is the landing point of the external ray $\mathcal { R } _ { - 2 } ( \theta )$ . For an arbitrary connected filled Julia set $K _ { P }$ , we have the following result of Carathéodory: the inverse $\psi _ { P }$ of the Böttcher map has a continuous extension from w : $| \boldsymbol { w } | > 1 \}$ to $\{ w : | w | \geqslant 1 \}$ if and only if $K _ { P }$ is locally connected. To understand what this means, imagine a set that is shaped like a comb. From any point in this set to any other point there is a continuous path that lies in the set, but it is possible for the two points to be very close and for the shortest path to be very long. This happens, for example, if the two points are the ends of neighboring teeth of the comb. A connected set X is called locally connected if every point has arbitrarily small connected neighborhoods. It is possible to build comb-like sets (with infinitely many teeth) that contain points for which all connected neighborhoods have to be large. The filled Julia sets in the examples in figures 2–5 are locally connected, but there are examples of filled Julia sets that are not locally connected. When $K _ { P }$ is locally connected, then all external rays land, and the landing point is a continuous function of the argument. Under these circumstances, we have a natural and useful parametrization of the Julia set $J _ { P }$ .

# 2.8 The Mandelbrot Set M

We shall now restrict our attention to quadratic polynomials of the form $Q _ { c } ,$ . These are parametrized by the complex number c, and in this context we shall refer to the complex plane as the parameter plane, or c-plane. We would like to understand the family of dynamical systems that arise when we iterate the polynomials $Q _ { c } .$ Our goal will be to do this by dividing the c-plane into regions that correspond to polynomials with qualitatively the same dynamics. These regions will be separated by their boundaries, which together form the so-called bifurcation set. This consists of “unstable” c-values: that is, values of c for which there are other values arbitrarily nearby that give rise to qualitatively different dynamical behavior. In other words, a parameter c belongs to the bifurcation set if a small perturbation of c can make an important difference to the dynamics.

Recall the dichotomy that we stated earlier: the Julia set $J _ { c }$ is connected if the critical point 0 belongs to the filled Julia set $K _ { c }$ and is totally disconnected if 0 belongs to the attracting basin $A _ { c } ( \infty )$ . This dichotomy motivates the following definition: the Mandelbrot set M consists of the c-values for which $J _ { c }$ is connected. That is,

$$
M = \{c \in \mathbb {C} \mid Q _ {c} ^ {k} (0) \nrightarrow \infty \text {   as   } k \to \infty \}.
$$

Since the Julia set represents the chaotic part of the dynamical system given by $\displaystyle Q _ { c } ,$ , the dynamical behavior is certainly qualitatively affected by whether c belongs to M or not. We have therefore made a start toward our goal, but the division of the plane into M and C  M is very coarse, and it does not obviously give us the complete understanding we are looking for.

The important set is in fact not M, but its boundary ∂M, which is illustrated in figure 7. Notice that this set has a number of “holes” (in fact, infinitely many). The Mandelbrot set itself is obtained by filling in all these holes. More precisely, the complement of ∂M consists of an infinite collection of connected components, of which one, the outside of the set, stretches off to infinity, while all the others are bounded. The “holes” are the bounded components.

This definition is similar to the definition of the Julia set of a polynomial. It is easy to define the filled Julia set, and the Julia set is then defined as its boundary. The Julia set provides a lot of structure in the dynamical plane, the z-plane. The Mandelbrot set is similarly easy to define, and its boundary provides a lot of structure in the c-plane. Remarkably, even though each Julia set concerns just one dynamical system, while the Mandelbrot set concerns an entire family of systems, there are close analogies between them, as will become clear.

![](images/1267d5d63f1d62acde06de7040a149d8b5f89b9db5053c43ffd356af1a1fbd2d.jpg)

<details>
<summary>natural_image</summary>

Fractal geometric pattern with symmetrical fractal-like shapes (no text or symbols)
</details>

Figure 7 The boundary ∂M of the Mandelbrot set.

Pioneering work on holomorphic dynamics in general and quadratic polynomials in particular was carried out in the early 1980s by Adrien Douady and John H. Hubbard. They introduced the name “Mandelbrot set” and proved several results about it. In particular, they defined a sort of Böttcher map, denoted by $\varPhi _ { M }$ , for the Mandelbrot set, which is a map from the complement of the Mandelbrot set to the complement of the closed unit disk.

The definition of $\varPhi _ { M }$ is actually quite simple: for each c let $\phi _ { M } ( { \boldsymbol { c } } )$ equal $\varphi _ { c } ( c )$ , where $\varphi _ { c }$ is the Böttcher map for the parameter c. However, Douady and Hubbard did more than merely define $\phi _ { M } \mathbf { : }$ they proved that it is a holomorphic bijection with holomorphic inverse.

Once we have $\varPhi _ { M }$ we can make further definitions, just as we did with the Böttcher map. For instance, we can define a potential G on the complement of the Mandelbrot set by setting $G ( c ) = g _ { c } ( c ) = \log | \varPhi _ { M } ( c ) |$ . An equipotential is then a level set of $\varPhi _ { M }$ (that is, a set of the form $\{ c \ \in \ \mathbb { C } \ : \ | \phi _ { M } ( c ) | \ = \ r \}$ for some $r > 1 )$ and the external ray of argument θ is the set $\{ c \in \mathbb { C } : \arg ( \phi _ { M } ( c ) ) = 2 \pi \theta \}$ (that is, the inverse image of a radial line $\mathcal { R } _ { 0 } ( \theta ) )$ . The latter is denoted by ${ \mathcal { R } } _ { M } ( \theta )$ and it is asymptotic to the radial line of argument θ. The rational external rays are known to land (see figure 8).

It follows from the above that as t approaches zero, the equipotential of potential t, together with its interior, gets closer and closer to M: that is, M is the intersection of all such sets. Hence, M is a connected, closed, bounded subset of the plane.

![](images/5344daecba11ad32e5df37a8dd377a8c50c37c2e5a8a753b80672d8556ebfc12.jpg)

<details>
<summary>natural_image</summary>

Symmetrical fractal pattern with concentric circles and radial lines, no text or symbols present
</details>

Figure 8 Some equipotentials of M and the external rays of arguments θ of periods 1, 2, 3, and 4. In counterclockwise direction the arguments between 0 and $\begin{array} { l } { { \frac { 1 } { 2 } } } \end{array}$ are $0 , { \frac { 1 } { 1 5 } } , { \frac { 2 } { 1 5 } } ,$ , 17 , 31 5 , 41 5 , 27 , 13 , 61 5 , 37 , a n d 71 ${ \begin{array} { l } { { \frac { 1 } { 7 } } , \ { \frac { 3 } { 1 5 } } , \ { \frac { 4 } { 1 5 } } , \ { \frac { 2 } { 7 } } , \ { \frac { 1 } { 3 } } , \ { \frac { 6 } { 1 5 } } , \ { \frac { 5 } { 7 } } } \end{array} }$ ${ \frac { 7 } { 1 5 } } ;$ 2 15 15 ; and symmetrically in clockwise direction they are $1 - \theta$ with θ as above. The external rays of argument $\frac { 1 } { 7 }$ and $\frac { 2 } { 7 }$ are landing at the root point of the hyperbolic component that has $c _ { 0 } ,$ , the parameter value of the Douargument $\frac { 3 } { 1 5 }$ y raband $\frac { 4 } { 1 5 }$ in figure 2, as its center. The rays ofare landing at the root point of the

# 2.8.1 J-Stability

As we have mentioned and as figure 7 suggests, the complement of ∂M has infinitely many connected components. These components are of great dynamical significance: if c and $c ^ { \prime }$ are two parameters taken from the same component, then the dynamical systems arising from $Q _ { c }$ and $Q _ { c ^ { \prime } }$ can be shown to be essentially the same. To be precise, they are J-equivalent, which means that there is a continuous change of variables that converts the dynamics on one Julia set to the dynamics on the other. If c belongs to the boundary ∂M, then there are parameter values $c ^ { \prime }$ arbitrarily close to c for which $Q _ { c }$ and $Q _ { c ^ { \prime } }$ are not J-equivalent, so ∂M is the “bifurcation set with respect to J-stability.” We shall comment on the global structural stability later.

# 2.8.2 Hyperbolic Components

From now on, we shall use the word “component” to refer to the holes of the Mandelbrot set—that is, to the bounded components of the complement of ∂M.

We start by considering the component containing $c \ = \ 0$ , the central component ${ \mathcal { H } } _ { 0 } .$ . Recall from section 2.3 that, after a suitable change of variables, one can change the polynomial $F _ { \lambda } ( z ) ~ = ~ \lambda z + z ^ { 2 }$ into the polynomial $\displaystyle Q _ { c } ,$ , where the parameters λ and c are related by the equation $\begin{array} { r } { c = \frac { 1 } { 2 } \lambda - \frac { 1 } { 4 } \lambda ^ { 2 } } \end{array}$ . The parameter λ has a dynamical meaning: the origin is a fixed point of $F _ { \lambda }$ and λ is its multiplier. This knowledge tells us that the corresponding $Q _ { c }$ has a fixed point of multiplier $\lambda ;$ we denote the fixed point by $\alpha _ { c } .$ . For $| \lambda | < 1$ the fixed point is attracting.

The unit disk $\{ \lambda : | \lambda | < 1 \}$ corresponds to the central component ${ \mathcal { H } } _ { 0 } ,$ , and the function that takes a parameter c in ${ \mathcal { H } } _ { 0 }$ to the corresponding parameter λ in the unit disk is called the multiplier map, and is denoted by $\rho _ { \mathcal { H } _ { 0 } }$ . Thus, $\rho _ { \mathcal { H } _ { 0 } } ( c )$ is the multiplier of the fixed point $\alpha _ { c }$ of the polynomial $Q _ { c }$ . The multiplier map $\rho \mathcal { H } _ { 0 }$ is a holomorphic isomorphism from ${ \mathcal { H } } _ { 0 }$ to the unit disk. As we have just seen, the inverse map is given by $\begin{array} { r } { \rho _ { \mathcal { H } _ { 0 } } ^ { - 1 } ( \lambda ) = \frac { 1 } { 2 } \lambda - \frac { 1 } { 4 } \lambda ^ { 2 } } \end{array}$ . This map extends continuously to the unit circle, and thereby gives us a parametrization of the boundary of the central component ${ \mathcal { H } } _ { 0 }$ by points λ of modulus 1. The image of the unit circle under the map $\lambda \mapsto { \textstyle { \frac { 1 } { 2 } } } \lambda - { \textstyle { \frac { 1 } { 4 } } } \lambda ^ { 2 }$ is a cardioid. This explains the heartlike shape of the largest part of the Mandelbrot set, which can be seen in figure 7.

Any quadratic polynomial has two fixed points if we count with multiplicity (in fact, two distinct ones unless $\begin{array} { r } { c = \frac { 1 } { 4 } ) } \end{array}$ . The central component ${ \mathcal { H } } _ { 0 }$ is characterized as the component of c-values for which $Q _ { c }$ has an attracting fixed point. For any c outside the cardioid, $Q _ { c }$ has two repelling fixed points, but it may have an attracting periodic orbit of a period greater than 1. It is an important fact that the attracting basin of an attracting periodic orbit always contains a critical orbit. Therefore, for any quadratic polynomial there can be at most one attracting periodic orbit.

We call a component  of the Mandelbrot set a hyperbolic component if, for every parameter c in $\mathcal { H }$ , the polynomial $Q _ { c }$ has an attracting periodic orbit. For any given hyperbolic component, the periods of the attracting periodic orbits will be the same. There is a corresponding multiplier map $\rho _ { \mathcal { H } }$ , from to the unit disk, which assigns to each parameter c in the multiplier of the attracting periodic orbit. This multiplier map is always a holomorphic isomorphism that extends continuously to the boundary ∂ of .

The points $\rho _ { \mathcal { H } } ^ { - 1 } ( 0 )$ and $\rho _ { \mathcal { H } } ^ { - 1 } ( 1 )$ are called the center and the root of . The center of is the unique c in for which the periodic orbit of $Q _ { c }$ is super-attracting. As for the root, if the period of the component is $k ,$ then it will be the landing point for a pair of external rays of periodic arguments of period k. (For the central component ${ \mathcal { H } } _ { 0 }$ there is only one ray assigned.) Conversely, every external ray with such an argument lands at the root point of a hyperbolic component of period k. Thus, the arguments of these rays give addresses to the hyperbolic components. This can be seen in figure 8, from which one can read off the mutual positions of all the components of periods 1–4.

As a consequence of the above, the number of hyperbolic components corresponding to a certain period k can be determined both as the number of roots in the polynomial $Q _ { c } ^ { k } ( 0 )$ that are not roots in $Q _ { c } ^ { \ell } ( 0 )$ ) for some $\ell \ < \ k$ and also as the number of pairs of rational arguments with denominator $2 ^ { k } - 1$ that cannot be expressed with denominator $2 ^ { \ell } - 1$ for some $\ell < k$ .

For any component  with center $c _ { 0 } \operatorname { l e t } \mathcal { R } _ { M } ( \theta _ { - } )$ ) and $\mathcal { R } _ { M } ( \theta _ { + } )$ be the pair of rays landing at the root point. Then, in the dynamical plane of $Q _ { c _ { 0 } } ,$ the pair of rays $\mathcal { R } _ { c _ { 0 } } ( \theta _ { - } )$ and $\mathcal { R } _ { c _ { 0 } } ( \theta _ { + } )$ are adjacent to the Fatou component of $Q _ { c _ { 0 } }$ containing $c _ { 0 }$ , and they land at the root point of that Fatou component.

# 2.8.3 Structural Stability

Suppose that $Q _ { c }$ has a super-attracting periodic orbit of period k, and let $z _ { 0 }$ be a point in this orbit. Then $Q _ { c } ^ { k } ( z _ { 0 } ) = z _ { 0 }$ , and the derivative of $Q _ { c } ^ { k }$ at $z _ { 0 }$ is 0. It follows from the chain rule that there is at least one $z _ { i }$ in the orbit at which the derivative of $Q _ { c }$ is 0: that is, 0 belongs to the orbit. Therefore, the center of a hyperbolic component cannot be structurally stable, since the critical orbit of the center-polynomial is finite, but it is infinite for all nearby polynomials. However, if we remove from the complex plane not just ∂M but also all the centers of hyperbolic components, then we obtain the splitting we have been looking for: any connected component of the remaining set forms a structurally stable region. For any pair of parameter values c and $c ^ { \prime }$ in such a component, $Q _ { c }$ and $Q _ { c }$ are conjugate, meaning that there is a continuous change of variables in the plane that converts the dynamics of one polynomial into those of the other.

# 2.8.4 Conjectures

The above discussion raises an obvious question: we have a good understanding of the hyperbolic components of the complement of ∂M, but are there components that are not hyperbolic? The following conjecture expresses a widely held belief, but it is as yet unproved.

The hyperbolicity conjecture. All the bounded components of the complement of ∂M are hyperbolic.

The hyperbolicity conjecture can be stated in greater generality for rational functions, where it says that every rational function can be approximated arbitrarily closely by a hyperbolic rational function. Here, “hyperbolic” means that the dynamics is expanding on the Julia set. We shall not go further into this, but only mention that the dynamics on the Julia set is expanding for every $Q _ { c }$ with c in a hyperbolic component of $M ,$ and also in the unbounded component, the complement of M. The Julia set $J _ { c }$ can in these cases be thought of as a “strange repeller”: the dynamics is chaotic and the geometry is fractal (except for $c = 0 )$ .

The main conjecture about the Mandelbrot set is, however, the following.

The local connectivity conjecture. The Mandelbrot set is locally connected.

This conjecture, often referred to as MLC, is important for many reasons. To begin with, it is known that it implies the hyperbolicity conjecture. Second, if M is locally connected, then $\Psi _ { M } ,$ , the inverse of $\varPhi _ { M }$ , which is a holomorphic bijection from the set outside the closed unit disk to the complement of the Mandelbrot set, has a continuous extension to the unit circle, and all external rays land in a continuous manner. This would give us a useful parametrization of ∂M. One can then give a beautifully simple abstract combinatorial description of $M ,$ despite the fact that ∂M is a complicated fractal. (Mitsuhiro Shishikura has proved that the hausdorff dimension [III.17] of ∂M is the maximum possible in the plane, namely 2.)

# 2.9 Universality of M

The Mandelbrot set is remarkably ubiquitous. For example, homeomorphic copies of M appear inside M itself, as is apparent from figure 9. Inside other families of holomorphic mappings that depend holomorphically on some parameter, we again find homeomorphic copies of M. For this reason, M is said to be universal. Douady and Hubbard have captured the reason behind the phenomenon of universality by defining a notion of a quadratic-like mapping. The kth iterate of a quadratic polynomial is globally a polynomial of degree $2 ^ { k }$ , but locally it may behave like a quadratic polynomial. The same is true for a rational function or an iterate of it. By a quadratic-like mapping we mean a triple $( f , V , W )$ where V and W are open simply connected domains (that ${ \mathrm { i } } s ,$ connected open sets without holes), $\bar { V } \subset W .$ , and $f$ is a holomorphic map that maps V onto W with degree 2. (This means that every point in W has two preimages, up to multiplicity, in V .) Such a map $f$ has a single critical point ω in $V ,$ , and behaves in many ways like a quadratic polynomial. The filled Julia set $K _ { f }$ is defined as the set of points z in V for which the iterates $f ^ { k } ( z )$ stay in V for all $k \geqslant 0 .$ A dichotomy similar to the one for quadratic polynomials holds for quadratic-like mappings as well: $K _ { f }$ is connected if and only if the critical point ω is contained in $K _ { f }$ . For any quadratic-like mapping with a connected filled Julia set, Douady and Hubbard have defined a strategy, called straightening, which associates with the mapping a unique c-value in M. For a family of quadratic-like mappings $\{ f _ { \lambda } \} _ { \lambda \in { \cal { N } } }$ the Mandelbrot set $M _ { \varLambda }$ is defined as the set of λ for which $K _ { f _ { \lambda } }$ is connected. We obtain through straightening a mapping $\Xi : M _ { \Lambda } \to M$ , which takes λ to the uniquely associated c-value.

![](images/5b4940fdbc88acf5618ca6a55daba1f2061ca8f784ad9572d6790a54bf56accc.jpg)

<details>
<summary>natural_image</summary>

Fractal pattern with branching structures and central void (no text or symbols)
</details>

Figure 9 A copy of M within M. The address of the copy is given by the arguments of the two external rays that+ land at the cusp, the root point of the copy. Here the arguments are $\frac { 3 } { 1 5 }$ and $\frac { 4 } { 1 5 }$ . Compare with figure 8. The rays are drawn to indicate where the “decorations” should be cut off in order to have the bare copy of M.

In the copy of M shown in figure 9, the “center” associated with $c = 0$ in M corresponds to a polynomial $Q _ { c _ { 0 } }$ for which the critical point 0 is periodic of period 4, and for which a suitable restriction of the fourth iterate $f _ { c _ { 0 } } = Q _ { c _ { 0 } } ^ { 4 }$ is quadratic-like from $V _ { 0 }$ to its image $W _ { 0 }$ . Moreover, there is a neighborhood $\gamma _ { 0 }$ of $c _ { 0 }$ in the c-plane such that for any $c$ in $\gamma _ { 0 }$ the restriction of $f _ { c } ~ = ~ Q _ { c } ^ { 4 }$ to V0 is a quadratic-like map from $V _ { 0 }$ to its image $W _ { c } ,$ , and such that the map Ξ is a homeomorphism from $M \gamma _ { 0 }$ to M .

The infinitely many copies of M that appear inside M may suggest that M has a self-similarity property. However, there is another phenomenon that pulls in the opposite direction. The c-values for which the critical point 0 is pre-periodic form a dense subset of ∂M. If c˜ is one of these special c-values, then there are two contexts in which one may look at magnifications of smaller and smaller neighborhoods of c˜: the first is the Julia set $J _ { \tilde { c } }$ of the polynomial $Q _ { \tilde { c } }$ in neighborhoods of $z = { \tilde { c } } ,$ and the second is the Mandelbrot set in neighborhoods of $c = { \tilde { c } } .$ . It turns out that the pictures are asymptotically similar, which means that the greater the magnification, and the smaller the neighborhood, the more similar the two pictures become.

This is an extraordinary fact. Indeed, it may even seem to be impossible, since in any neighborhood of c˜ the Mandelbrot set contains infinitely many copies of itself, while the Julia set is known to contain no such copies. The explanation for the apparent paradox is that the copies of the Mandelbrot set get smaller very quickly as their distance to c˜ decreases. Hence, if one magnifies a small enough neighborhood, the copies that are there are practically invisible.

# 2.10 Newton’s Method Revisited

Let us return briefly to Newton’s method for polynomials. Consider any polynomial P of degree d $\geqslant 2$ that has only simple roots. Then the Newton function $N _ { P }$ is a rational function of degree d, and each simple root of $P$ is a super-attracting fixed point of $N _ { P }$ . For quadratic polynomials the number of roots of P coincides with the number of critical points of $N _ { P }$ (since $2 d - 2 = 2$ when $d = 2 )$ . For polynomials of degree d $> 2$ there are more critical points than the roots can account for.

Cayley considered Newton’s method for quadratic polynomials with two distinct roots $P ( z ) = ( z - r _ { 1 } ) ( z -$ $r _ { 2 } )$ . He showed that the function $\mu ( z ) = ( z - r _ { 1 } ) / ( z$ $r _ { 2 } )$ , which maps the root r1 onto 0 and the root r2 onto , provides a change of variables that turns $N _ { P }$ into the quadratic polynomial $Q _ { 0 }$ on the Riemann sphere ${ \hat { \mathbb { C } } } .$ When one translates the dynamics of $Q _ { 0 }$ to the dynamics of Newton’s method one finds that the unit circle corresponds to the bisector of $r _ { 1 }$ and $r _ { 2 }$ and that all points in the half-plane containing $r _ { i } , i = 1 , 2 ,$ , are therefore attracted to $r _ { i }$ under iteration by $N _ { P }$ .

Cayley announced that he would write about Newton’s iteration for cubic polynomials. However, it took about a hundred years before any such paper appeared. For a cubic polynomial P with three simple roots the Newton function $N _ { P }$ has three super-attracting fixed points, each of which gives rise to an attracting basin. The Julia set of $N _ { P }$ is the common boundary of these three basins, and is therefore a complicated fractal set. Moreover, $N _ { P }$ has an extra critical point since $2 d - 2 = 4$ for $d = 3$ . The extra critical point may be attracted to one of the roots under iteration, or it can have its own independent behavior. In order to catch the behavior of all cubic polynomials under Newton’s iteration (except the one with one root of multiplicity three) it is sufficient to consider the one-parameter family of polynomials $P _ { \lambda } ( z ) = ( z - 1 ) ( z - { \textstyle { \frac { 1 } { 2 } } } - \lambda ) ( z - { \textstyle { \frac { 1 } { 2 } } } + \lambda )$ . The extra critical point for the corresponding Newton function $N _ { \lambda }$ then turns out to be at the origin. Suppose that we associate three colors, for instance red, blue, and green, with the three roots $1 , { \frac { 1 } { 2 } } + \lambda , { \frac { 1 } { 2 } }$ −λ. We can then color the λ-plane, which is the parameter plane in this context, as follows. A parameter value λ is colored red, blue, or green if the critical point 0 is attracted under iteration by $N _ { \lambda }$ to the root of that color. If it is not attracted to any of the three roots, then we color with a fourth color, yellow, say. The universality of the Mandelbrot set is thereby demonstrated: in the λ-plane one can observe yellow copies of ${ \mathrm { i t } } ,$ which one can explain by showing that families of suitably restricted iterates of $N _ { \lambda }$ are quadratic-like.

# 3 Concluding Remarks

We have illustrated several results in holomorphic dynamics through examples, including the transferring of definitions and results from the dynamical planes to the parameter plane. The structures of the filled Julia sets and the Mandelbrot set are partly understood through analysis of their complements, linked together via the Böttcher maps $\varphi _ { c }$ and $\varPhi _ { M }$ . The functions that are used for changing variables in J-stability and structural stability are examples of so-called quasiconformal mappings. This is a concept that was introduced into holomorphic dynamics in the early 1980s by

Dennis Sullivan. They are indispensable for discussing change of complex structure, straightening, holomorphic motion, surgery, and many other phenomena. The interested reader is referred to the books listed below. The first two contain expository papers, the third is a graduate textbook, and the fourth is a collection of papers. They all contain many further references.

Acknowledgments. The computer drawings in this article were obtained from a program written by Christian Henriksen.

# Further Reading

Devaney, R. L., and L. Keen, eds. 1989. Chaos and Fractals. The Mathematics Behind the Computer Graphics. Proceedings of Symposia in Applied Mathematics, volume 39. Providence, RI: American Mathematical Society.   
. 1994. Complex Dynamical Systems. The Mathematics Behind the Mandelbrot and Julia Sets. Proceedings of Symposia in Applied Mathematics, volume 49. Providence, RI: American Mathematical Society.   
Lei, T., ed. 2000. The Mandelbrot Set, Theme and Variations. London Mathematical Society Lecture Note Series, volume 274. Cambridge: Cambridge University Press.   
Milnor, J. 1999. Dynamics in One Complex Variable. Weisbaden: Vieweg.

# IV.15 Operator Algebras

# Nigel Higson and John Roe

# 1 The Beginnings of Operator Theory

We can ask two basic questions about any equation, or system of equations: is there a solution, and, if there is, is it unique? Experience with finite systems of linear equations indicates that the two questions are interconnected. Consider for instance the equations

$$
2 x + 3 y - 5 z = a,
$$

$$
x - 2 y + z = b,
$$

$$
3 x + y - 4 z = c.
$$

Notice that the left-hand side of the third equation is the sum of the left-hand sides of the first two. As a result, no solution to the system exists unless $a + b = c .$ . But if $a + b \ = \ c ,$ then any solution of the first two equations is also a solution of the third; and in any linear system involving more unknowns than equations, solutions, when they exist, are never unique. In the present case, if $( x , y , z )$ is a solution, then so is $( x + t , y + t , z + t )$ , for any t. Thus the same phenomenon (a linear relation among the equations) that prevents the system from admitting solutions in some cases also prevents solutions from being unique in other cases.

To make the relation between existence and uniqueness of solutions more precise, consider a general system of linear equations of the form

$$
k _ {1 1} u _ {1} + k _ {1 2} u _ {2} + \dots + k _ {1 n} u _ {n} = f _ {1},
$$

$$
k _ {2 1} u _ {1} + k _ {2 2} u _ {2} + \dots + k _ {2 n} u _ {n} = f _ {2},
$$

$$
k _ {n 1} u _ {1} + k _ {n 2} u _ {2} + \dots + k _ {n n} u _ {n} = f _ {n}
$$

consisting of n equations in n unknowns. The scalars $k _ { j i }$ form a matrix of coefficients and the problem is to solve for the $u _ { i }$ in terms of the $f _ { j }$ . The general theorem illustrated by our particular numerical example above is that the number of linear conditions that the $f _ { j }$ must satisfy if a solution is to exist is equal to the number of arbitrary constants appearing in the general solution when a solution does exist. To use a more technical vocabulary, the dimension of the kernel [I.3 §4.1] of the matrix $K \ = \ \{ k _ { j i } \}$ is equal to the dimension of its cokernel. In the example, these numbers are both 1.

A little more than a hundred years ago, fredholm [VI.66] made a study of integral equations of the type

$$
u (y) - \int k (y, x) u (x) \mathrm{d} x = f (y).
$$

These arose from questions in theoretical physics, and the problem was to solve for the function u in terms of the function $f .$ Since an integral can be thought of as a limit of finite sums, Fredholm’s equation is an infinitedimensional counterpart of the finite-dimensional linear systems considered above, in which vectors with n components are replaced by functions with values at infinitely many different points x. (Strictly speaking, Fredholm’s equation is analogous to a matrix equation of the type $u - K u = f$ rather than $K u = f .$ The altered form of the left-hand side has no effect on the overall behavior of the matrix equation, but it does considerably alter the behavior of the integral equation. As we shall see, Fredholm was fortunate to work with a class of equations whose behavior mirrors that of matrix equations very closely.)

A very simple example is

$$
u (y) - \int_ {0} ^ {1} u (x) \mathrm{d} x = f (y).
$$

To solve this equation, it helps to observe that the quantity $\textstyle \int _ { 0 } ^ { 1 } u ( x ) \mathrm { d } x$ , when thought of as a function of $_ { \mathcal { V } , }$ is a constant. Thus in the homogeneous case $( f \equiv 0 )$ , the only possible solutions for u(y) are the constant functions. On the other hand, for a general function $f ,$ solutions exist if and only if the single linear condition $\int _ { 0 } ^ { 1 } f ( y ) { \mathrm { d } } y = 0$ is satisfied. So in this example the dimension of the kernel and the dimension of the cokernel are both 1. Fredholm set out on a systematic exploration of the analogy between matrix theory and integral equations that this example suggests. He was able to prove that, for equations of his type, the dimensions of the kernel and of the cokernel are always finite and equal.

Fredholm’s work sparked the imagination of hilbert [VI.63], who made a detailed study of the integral operators that transform $u ( y )$ into $\int k ( y , x ) u ( x ) \mathrm { d } x$ , in the special case where the real-valued function k is symmetric, meaning that $k ( x , y ) ~ = ~ k ( y , x )$ . The finite-dimensional counterpart of Hilbert’s theory is the theory of real symmetric matrices. Now if K is such a matrix, then a standard result from linear algebra asserts that there is an orthonormal basis consisting of eigenvectors [I.3 §4.3] for K, or equivalently that there is a unitary matrix U such that $U ^ { - 1 } T U$ is diagonal. (Unitary means that U is invertible and preserves the lengths of vectors: $\| U \nu \| = \| \nu \|$ for all vectors v.) Hilbert obtained an analogous theory for all symmetric integral operators. He showed that there exist functions $u _ { 1 } ( y ) , u _ { 2 } ( y ) , . .$ . and real numbers $\lambda _ { 1 } , \lambda _ { 2 } , \ldots$ such that

$$
\int k (y, x) u _ {n} (x) \mathrm{d} x = \lambda_ {n} u _ {n} (y).
$$

Thus $u _ { n } ( y )$ is an eigenfunction for the integral operator, with eigenvalue $\lambda _ { n }$ .

In most cases it is hard to calculate $u _ { n }$ and $\lambda _ { n }$ explicitly, but calculation is possible when $k ( x , y ) =$ $\phi ( x - y )$ for some periodic function φ. If the range of integration is [0, 1] and the period of φ is 1, then the eigenfunctions are cos(2kπy), $k = 0 , 1 , 2 , \ldots$ , and sin(2kπy), $k \ = \ 1 , 2 , \ldots$ In this case, the theory of fourier series [III.27] tells us that a general function $f ( y )$ on [0, 1] can be expanded as the sum of a series $\begin{array} { r } { \sum ( a _ { k } \cos 2 k \pi y + b _ { k } \sin 2 k \pi y ) } \end{array}$ of cosines and sines. Hilbert showed that, in general, there is an analogous expansion

$$
f (y) = \sum a _ {n} u _ {n} (y)
$$

in terms of the eigenfunctions for any symmetric integral operator. In other words, the eigenfunctions form a basis, just as in the finite-dimensional case. Hilbert’s result is now called the spectral theorem for symmetric integral operators.

# 1.1 From Integral Equations to Functional Analysis

Hilbert’s spectral theorem led to an explosion of activity, since integral operators arise in many different areas of mathematics (including, for example, the dirichlet problem [IV.12 §1] in partial differential equations and the representation theory of compact groups [IV.9 §3]). It was soon recognized that these operators are best viewed as linear transformations on the hilbert space [III.37] of all functions $u ( y )$ such that $\int | u ( y ) | ^ { 2 } \mathrm { d } y < \infty$ . Such functions are called square-integrable, and the collection of all of them is denoted $L ^ { 2 } [ 0 , 1 ]$ .

With the important concept of Hilbert space available, it became convenient to examine a much broader range of operators than the integral operators initially considered by Fredholm and Hilbert. Since Hilbert spaces are vector spaces [I.3 §2.3] and metric spaces [III.56], it made sense to look first at operators from a Hilbert space to itself that are both linear and continuous: these are usually called bounded linear operators. The analogue of the symmetry condition $k ( x , y ) ~ =$ $k ( y , x )$ on integral operators is the condition that a bounded linear operator T be self-adjoint, which is to say that $\langle T u , \nu \rangle = \langle u , T \nu \rangle$ for all vectors u and v in the Hilbert space (the angle brackets denote the inner product). A simple example of a self-adjoint operator is the multiplication operator by a real-valued function $m ( \boldsymbol { y } ) ;$ this is the operator M defined by the formula $( M u ) ( y ) = m ( y ) u ( y )$ . (The finite-dimensional counterpart to a multiplication operator is a diagonal matrix $K ,$ which multiplies the jth component of the vector by the matrix entry $k _ { j j } . )$

Hilbert’s spectral theorem for symmetric integral operators tells us that every such operator can be given a particularly nice form: with respect to a suitable $\mathrm { \ " b a s i s \ " }$ of $L ^ { 2 } [ 0 , 1 ]$ , namely a basis of eigenfunctions, it will have an infinite diagonal matrix. Moreover, the basis vectors can be chosen to be orthogonal to each other. For a general self-adjoint operator, this is not true. Consider, for instance, the multiplication operator from $L ^ { 2 } [ 0 , 1 ]$ to itself that takes each square-integrable function $u ( y )$ to the function yu(y). This operator has no eigenvectors [I.3 §4.3], since if λ is an eigenvalue [I.3 §4.3], then we need $y u ( y ) ~ = ~ \lambda u ( y )$ for every y, which implies that $u ( y ) = 0$ for every y not equal to λ, and hence that $\begin{array} { r } { \int | u ( y ) | ^ { 2 } \mathrm { d } y = 0 } \end{array}$ . However, this example is not particularly worrying, since a multiplication operator of this kind is a sort of continuous analogue of the operator defined by a diagonal matrix.

It turns out that if we enlarge our concept of “diagonal” to include multiplication operators, then all self-adjoint operators are “diagonalizable,” in the sense that, after a suitable “change of basis,” they become multiplication operators.

To make this statement precise, we need the notion of the spectrum [III.86] of an operator T . This is the set of complex numbers λ for which the operator $T - \lambda I$ does not have a bounded inverse (here I is the identity operator on Hilbert space). In finite dimensions the spectrum is precisely the set of eigenvalues, but in infinite dimensions this is not always so. Indeed, whereas every symmetric matrix has at least one eigenvalue, a self-adjoint operator, as we have just seen, need not. As a result of this, the spectral theorem for bounded self-adjoint operators is phrased not in terms of eigenvalues but in terms of the spectrum. One way of formulating it is to state that any self-adjoint operator T is unitarily equivalent to a multiplication operator $( M u ) ( y ) = m ( y ) u ( y )$ , where the closure of the range of the function m(y) is the spectrum of T . Just as in the finite-dimensional case, a unitary operator is an invertible operator U that preserves the lengths of vectors. To say that T and M are unitarily equivalent is to say that there is some unitary map $U ,$ which we can think of as an analogue of a change-of-basis matrix, such that T $U ^ { - 1 } M U$ . This generalizes the statement that any real symmetric matrix is unitarily equivalent to a diagonal matrix with the eigenvalues along the diagonal.

# 1.2 The Mean Ergodic Theorem

A beautiful application of the spectral theorem was found by von neumann [VI.91]. Imagine a checkerboard on which are distributed a certain number of checkers. Imagine that for each square there is designated a “successor” square (in such a way that no two squares have the same successor), and that every minute the checkers are rearranged by moving each one to its successor square. Now focus attention on a single square and each minute record with a 1 or 0 whether or not there is a piece on the square. This produces a succession of readings $R _ { 1 } , R _ { 2 } , R _ { 3 } , \ldots$ . like this:

$$
0 0 1 0 0 1 1 0 0 1 0 1 1 0 1 0 0 1 0 0 \dots .
$$

We might expect that over time, the average number of positive readings $R _ { j } = 1$ will converge to the number of pieces on the board divided by the number of squares. If the rearrangement rule is not complicated enough, then this will not happen. For example, in the most extreme case, if the rule designates each square as its own successor, then the readout will be either 00000  or 111111  , depending on whether or not we chose a square with a piece on it to begin with. But if the rule is sufficiently complicated, then the “time aver-$\begin{array} { r } { \mathrm { a g e } ^ { \prime \prime } \left( 1 / n \right) \sum _ { j = 1 } ^ { n } R _ { j } } \end{array}$ will indeed converge to the number of pieces on the board divided by the number of squares, as expected.

The checkerboard example is elementary, since in fact the only “sufficiently complicated” rules in this finite case are cyclic permutations of the squares of the board, and thus all the squares move past our observation post in succession. However, there are related examples where one observes only a small fraction of the data. For instance, replace the set of squares on a checkerboard with the set of points on a circle, and in place of the checkers, imagine that a subset S of a circle is marked as occupied. Let the rearrangement rule be the rotation of points on the circle through some irrational number of degrees. Stationed at a point x of the circle, we record whether x belongs to S, the first rotated copy of S, the second rotated copy of S, and so on to obtain a sequence of 0 or 1 readings as before. One can show that (for nearly every x) the time average of our observations will converge to the proportion of the circle occupied by S.

Similar questions about the relationship between time and space averages had arisen in thermodynamics and elsewhere, and the expectation that time and space averages should agree when the rearrangement rule is sufficiently complex became known as the ergodic hypothesis.

Von Neumann brought operator theory to bear on this question in the following way. Let H be the Hilbert space of functions on the squares of the checkerboard, or the Hilbert space of square-integrable functions on the circle. The rearrangement rule gives rise to a unitary operator U on H by means of the formula

$$
(U f) (y) = f (\phi^ {- 1} (y)),
$$

where φ is the function describing the rearrangement. Von Neumann’s ergodic theorem asserts that if no nonconstant function in H is fixed by U (this is one way of saying that the rearrangement rule is “sufficiently complicated”), then, for every function $f \in H$ , the limit

$$
\lim _ {n \rightarrow \infty} \frac {1}{n} \sum_ {j = 1} ^ {n} U ^ {j} f
$$

exists and is equal to the constant function whose value everywhere is the average value of f . (To apply this to our examples, take $f ( x )$ to be the function that is 1 if the point x is occupied and 0 otherwise.)

Von Neumann’s theorem can be deduced from a spectral theorem for unitary operators that is analogous to the spectral theorem for self-adjoint operators. Every unitary operator can be reduced to a multiplication operator, not by real-valued functions but by functions whose values are complex numbers of absolute value 1. The key to the proof then becomes a statement about complex numbers of absolute value 1: if z is such a complex number, different from 1, then the expression $\textstyle ( 1 / n ) \sum _ { j = 1 } ^ { n } z ^ { j }$ approaches zero as n → ∞. This in turn is easily proved using the formula for the sum of a geometric series, ${ \textstyle \sum _ { j = 1 } ^ { n } z ^ { j } = z ( 1 - z ^ { n } ) / ( 1 - z ) }$ . (More detail can be found in ergodic theorems [V.9].)

# 1.3 Operators and Quantum Theory

Von Neumann realized that Hilbert spaces and their operators provide the correct mathematical tools to formalize the laws of quantum mechanics, introduced in the 1920s by Heisenberg and Schrödinger.

The state of a physical system at any given instant is the list of all the information needed to determine its future behavior. If, for instance, the system consists of a finite number of particles, then classically its state consists of the list of the position and momentum vectors of all the constituent particles. By contrast, in von Neumann’s formulation of quantum mechanics one associates with each physical system a Hilbert space H, and a state of the system is represented by a unit vector u in H. (If u and v are unit vectors and v is a scalar multiple of u, then u and v determine the same state.)

Associated with each observable quantity (perhaps the total energy of the system, or the momentum of one particle within the system) is a self-adjoint operator Q on H whose spectrum is the set of all observed values of that quantity (hence the origin of the term “spectrum”). States and observables are related as follows: when a system is in the state described by a unit vector $u \in H ,$ , the expected value of the observable quantity corresponding to a given self-adjoint operator Q is the inner product Qu, u . This may not be a value that is ever actually measured: rather, it is the average of values that are obtained from many repeated experiments with the system when it is in the given state u. The relation between states and observables reflects the paradoxical behavior of quantum mechanics: it is possible, and in fact typical, for a system to exist in a “superposed” state, under which repeated identical experiments produce distinct outcomes. A measurement of an observable quantity will produce a determinate outcome if and only if the state of the system is an eigenvector for the operator associated with that quantity.

A distinctive feature of quantum theory is that the operators associated with different observables typically do not commute with one another. If two operators do not commute, then they will typically have no eigenvectors in common, and, as a result, simultaneous measurements of two different observables will typically not result in determinate values for both of them. A famous example is provided by the operators P and Q associated with the position and momentum of a particle moving along a line. They satisfy the Heisenberg commutation relation

$$
Q P - P Q = \mathrm{i} \hbar I,
$$

where  is a certain physical constant. (This is an instance of a general principle which relates the noncommutativity of observables in quantum mechanics to the Poisson bracket of the corresponding observables in classical mechanics: see mirror symmetry [IV.16 §§2.1.3, 2.2.1].) As a result, it is impossible for the particle simultaneously to have a determinate momentum and position. This is the uncertainty principle.

It turns out that there is an essentially unique way of representing the Heisenberg commutation relation using self-adjoint operators on Hilbert space: the Hilbert space H must be $L ^ { 2 } ( \mathbb { R } ) ;$ ; the operator P must be i d/dx; and the operator Q must be multiplication by x. This theorem allows one to determine explicitly the observable operators for simple physical systems. For example, in a system consisting of a particle on a line subject to a force directed toward the origin which is proportional to the distance from the origin (as if the particle were attached to a spring, anchored at the origin), the operator for total energy is

$$
E = - \frac {\hbar^ {2}}{2 m} \frac {\mathrm{d} ^ {2}}{\mathrm{d} x ^ {2}} + \frac {k}{2} x ^ {2},
$$

where k is a constant which determines the overall strength of the force. The spectrum of this operator is the set

$$
\{(n + \frac {1}{2}) \hbar (k / m) ^ {1 / 2}: n = 0, 1, 2, \dots \}.
$$

These are therefore the possible values for the total energy of the system. Notice that the energy can assume only a discrete set of values. This is another characteristic and fundamental feature of quantum theory.

Another important example is the operator of total energy for the hydrogen atom. Like the operator above, this may be realized as a certain explicit partial differential operator. It can be shown that the eigenvalues of this operator form a sequence proportional to $\{ - 1 , - { \frac { 1 } { 4 } } , - { \frac { 1 } { 9 } } , \dots \}$ . A hydrogen atom, when disturbed, may release a photon, resulting in a drop in its total energy. The released photon will have energy equal to the difference between the energies of the initial and final states of the atom, and therefore it is proportional to a number of the form $1 / n ^ { 2 } - 1 / m ^ { 2 }$ . When light from hydrogen is passed through a prism or diffraction grating, bright lines are indeed observed at wavelengths corresponding to these possible energies. Spectral observations of this sort provide experimental confirmation for quantum mechanical predictions.

So far we have discussed states of a quantum system only at a single instant. However, quantum systems evolve in time, just as classical systems do: to describe this evolution we need a law of motion. The time evolution of a quantum system is represented by a family of unitary operators $U _ { t } : H \to H ,$ , parametrized by the real numbers. If the system is in an initial state $u ,$ it will be in the state $U _ { t } u$ after t units of time. Because the passage of s units of time followed by t further units is the same as the passage of $s + t$ units, the unitary operators $U _ { t }$ satisfy the group law $U _ { s } U _ { t } = U _ { s + t } .$ An important theorem of Marshall Stone asserts that there is a one-to-one correspondence between unitary groups $\{ U _ { t } \}$ and self-adjoint operators E given by the formula

$$
\mathrm{i} E = \left(\frac {\mathrm{d} U _ {t}}{\mathrm{d} t}\right) _ {t = 0} = \lim _ {t \rightarrow 0} \frac {1}{t} (U _ {t} - I).
$$

The quantum law of motion is that the generator E corresponding in this way to time evolution is the operator associated with the observable “total energy.” When E is realized as a differential operator on a Hilbert space of functions (as in the examples above), this statement becomes a differential equation, the Schrödinger equation.

# 1.4 The GNS Construction

The time-evolution operators $U _ { t }$ of quantum mechanics satisfy the law $U _ { s } U _ { t } = U _ { s + t }$ . More generally, we define a unitary representation of a group [I.3 §2.1] G to be a family of unitary operators $U _ { g } ,$ one for each $g \in G ,$ satisfying the law $U _ { g _ { 1 } g _ { 2 } } = U _ { g _ { 1 } } U _ { g _ { 2 } }$ for all $g _ { 1 } , g _ { 2 } \in G . \mathrm { O r i g i - }$ nally introduced by frobenius [VI.58] as a tool for the study of finite groups, representation theory [IV.9] has become indispensable in mathematics and physics wherever the symmetries of a system must be taken into account.

If U is a unitary representation of G and v is a vector, then $\sigma : g \mapsto \langle U _ { g } \nu , \nu \rangle$ is a function defined on G. The law $U _ { g _ { 1 } g _ { 2 } } = U _ { g _ { 1 } } U _ { g _ { 2 } }$ implies that $\sigma$ has an important positivity property, namely

$$
\sum_ {g _ {1}, g _ {2} \in G} \overline {{a _ {g _ {1}}}} a _ {g _ {2}} \sigma (g _ {1} ^ {- 1} g _ {2}) = \left| \left| \sum a _ {g} U _ {g} v \right| \right| ^ {2} \geqslant 0,
$$

for any scalars $a _ { g } \in \mathbb { C } . { \mathrm { ~ A ~ } }$ function defined on G and having this positivity property is said to be positive definite. Conversely, from a positive-definite function one can build a unitary representation. This GNS construction (in honor of Israel Gelfand, Mark Naimark, and Irving Segal) begins by considering the group elements themselves as basis vectors in an abstract vector space. We can attempt to define an inner product on this vector space by means of the formula

$$
\langle g _ {1}, g _ {2} \rangle = \sigma (g _ {1} ^ {- 1} g _ {2}).
$$

The resulting object may differ from a genuine Hilbert space in two respects. First, there may be nonzero vectors whose length, as measured by the inner product, is zero (although the hypothesis that σ is positive definite does rule out the possibility that there might be vectors of negative length). Second, the completeness axiom [III.62] of Hilbert space theory may not be satisfied. However, there is a “completion” procedure which fixes both these deficiencies. Applied in the present case, it produces a Hilbert space $H _ { \sigma }$ that carries a unitary representation of G.

Versions of the GNS construction arise in several areas of mathematics. They have the advantage that the functions on which the constructions are based are easy to manipulate. For instance, convex combinations of positive-definite functions are again positive definite, and this allows geometrical methods to be applied to the study of representations.

# 1.5 Determinants and Traces

The original works of Fredholm and Hilbert borrowed heavily from traditional concepts of linear algebra, and in particular the theory of determinants [III.15]. In view of the complicated definition of the determinant even for finite matrices, it is perhaps not surprising that the infinite-dimensional situation presented extraordinary challenges. Very soon, much simpler alternative approaches were found that avoided determinants altogether. But it is interesting to note that the determinant, or to be more exact the related notion of the trace, has played an important role in recent developments on which we will report later in this article.

The trace of an $n \times n$ matrix is the sum of its diagonal entries. As with the determinant, the trace of a matrix A is equal to the trace of $B A B ^ { - 1 }$ for any invertible matrix B. In fact, the trace is related to the determinant by the formula det(exp $( A ) ) = \exp ( \operatorname { t r } ( A ) )$ (because of the invariance properties of trace and determinant, it is enough to check this for diagonal matrices, where it is easy). In infinite dimensions the trace need not make sense since the sum of the diagonal entries of an $\infty \times \infty$ matrix may not converge. (The trace of the identity operator is a case in point: the diagonal entries are all 1, and if there are infinitely many of them, then their sum is not well-defined.) One way to address this problem is to limit oneself to operators for which the sum is well-defined. An operator T is said to be of trace class if, for every two sequences $\{ u _ { j } \}$ and $\{ \nu _ { j } \}$ of pairwise orthogonal vectors of length 1, the sum $\textstyle \sum _ { j = 1 } ^ { \infty } \langle T u _ { j } , \nu _ { j } \rangle$ is absolutely convergent. A trace-class operator T has a well-defined and finite trace, namely the sum $\textstyle \sum _ { j = 1 } ^ { \infty } \langle T u _ { j } , u _ { j } \rangle$ (which is independent of the choice of orthonormal basis $\{ u _ { j } \} )$ .

Integral operators such as those appearing in Fredholm’s equation provide natural examples of traceclass operators. If $k ( y , x )$ is a smooth function, then the operator $\begin{array} { r } { T u ( y ) = \int k ( y , x ) u ( x ) } \end{array}$ dx is of trace class, and its trace is equal to $\textstyle \int k ( x , x )$ dx, which can be regarded as the “sum” of the diagonal elements of the “continuous matrix” k.

# 2 Von Neumann Algebras

The commutant of a set S of bounded linear operators on a Hilbert space H is the collection $S ^ { \prime }$ of all operators on H that commute with every operator in the set S. The commutant of any set is an algebra of operators on H. That is, if T1 and $T _ { 2 }$ are in the commutant, then so are $T _ { 1 } T _ { 2 }$ and any linear combination $a _ { 1 } T _ { 1 } + a _ { 2 } T _ { 2 }$ .

As mentioned in the previous section, a unitary representation of a group G on a Hilbert space H is a collection of unitary operators $U _ { g } ,$ labeled by elements of G, with the property that for any two group elements $g _ { 1 }$ and $g _ { 2 }$ the composition $U _ { g _ { 1 } } U _ { g _ { 2 } }$ is equal to $U _ { g _ { 1 } g _ { 2 } }$ . A von Neumann algebra is any algebra of operators on a complex Hilbert space H which is the commutant of some unitary representation of a group on H. Every von Neumann algebra is closed under adjoints and under limits of nearly every sort. For example, it is closed under pointwise limits: if $\left\{ T _ { n } \right\}$ is a sequence of operators in a von Neumann algebra M, and if $T _ { n } \nu $ T v, for every vector $\nu \in H$ , then $T \in M$ .

It is easy to check that every von Neumann algebra M is equal to its own double commutant $M ^ { \prime \prime }$ (the commutant of the commutant of M). Von Neumann proved that if a self-adjoint algebra M of operators is closed under pointwise limits, then M is equal to the commutant of the group of unitary operators in its commutant, and is therefore a von Neumann algebra.

# 2.1 Decomposing Representations

Let $g  U _ { g }$ be a unitary representation of a group G on a Hilbert space H. If a closed subspace H0 of H is mapped into itself by all the operators $U _ { g }$ , then it is said to be an invariant subspace for the representation. If $H _ { 0 }$ is invariant, then since the operators $U _ { g }$ map $H _ { 0 }$ to itself, their restrictions to $H _ { 0 }$ constitute another representation of G, called a subrepresentation of the original.

A subspace $H _ { 0 }$ is invariant for a representation, and so determines a subrepresentation, if and only if the orthogonal projection operator $P : H \to H _ { 0 }$ belongs to the commutant of that representation. This points to a close connection between subrepresentations and von Neumann algebras. In fact, von Neumann algebra theory can be thought of as the study of the ways in which unitary representations can be decomposed into subrepresentations.

A representation is irreducible if it has no nontrivial invariant subspace. A representation that does have a nontrivial invariant subspace $H _ { 0 }$ can be divided into two subrepresentations: those associated with $H _ { 0 }$ and those associated with its orthogonal complement $H _ { 0 } ^ { \perp }$ . Unless both the representations $H _ { 0 }$ and $H _ { 0 } ^ { \perp }$ are irreducible, we will be able to divide one or both of them into still smaller pieces by repeating the process that was just carried out for H. If the initial Hilbert space H is finite dimensional, then continuing in this way we will eventually decompose it into irreducible subrepresentations. In the language of matrices, we will obtain a basis for H with respect to which all the operators in the group are simultaneously block diagonal, in such a way that each block represents an irreducible group of unitary operators on a smaller Hilbert space.

Reducing a unitary representation on a finite-dimensional Hilbert space into irreducible subrepresentations is a bit like decomposing an integer into a product of prime factors. As with prime factorization, the decomposition process for a finite-dimensional unitary representation has only one possible end: there is, up to ordering, a unique list of irreducible representations into which a given unitary representation decomposes.

But in infinite dimensions the decomposition process faces a number of difficulties, the most surprising of which is that there may be two decompositions of the same representation into entirely different sets of irreducible subrepresentations.

In the face of this, a different form of decomposition suggests itself, which is roughly analogous to the factorization of an integer into prime powers instead of individual primes. Let us refer to the prime powers into which an integer is decomposed as its components. They have two characteristic properties: no two components share a common factor, and any two (proper) factors of the same component do share a common factor. Similarly, one can decompose a unitary representation into isotypical components, which have analogous properties: no two distinct isotypical components share a common (meaning isomorphic) subrepresentation, and any two subrepresentations of the same isotypical component have themselves a common subsubrepresentation. Any unitary representation (finite dimensional or not) can be decomposed into isotypical components, and this decomposition is unique.

In finite dimensions, every isotypical representation decomposes into a (finite) number of identical irreducible subrepresentations (like the prime factors of a prime power). In infinite dimensions this is not so. In effect, much of von Neumann algebra theory is concerned with analyzing the many possibilities that arise.

# 2.2 Factors

The commutant of an isotypical unitary representation is called a factor. Concretely, a factor is a von Neumann algebra M whose center, the set of all operators in M that commute with every member of M, consists of nothing more than scalar multiples of the identity operator. This is because projections in the center of M correspond to projections onto combinations of isotypical subrepresentations. Every von Neumann algebra can be uniquely decomposed into factors.

A factor is said to be of type I if it arises as the commutant of an isotypical representation that is a multiple of a single irreducible representation. Every type I factor is isomorphic to the algebra of all bounded operators on a Hilbert space. In finite dimensions, every factor is of type I, since as we already noted every isotypical representation decomposes into a multiple of one irreducible representation.

The existence of unitary representations with more than one decomposition into irreducible components is related to the existence of factors that are not of type I. Von Neumann, together with Francis Murray, investigated this possibility in a series of papers that mark the foundation of operator algebra theory. They introduced an order structure on the collection of subrepresentations of a given isotypical representation or, to put it in terms of the commutant, on the collection of projections in a given factor. If $H _ { 0 }$ and $H _ { 1 }$ are subrepresentations of the isotypical representation $H ,$ then we write $H _ { 0 } ~ \preceq ~ H _ { 1 }$ if $H _ { 0 }$ is isomorphic to a subrepresentation of $H _ { 1 }$ . Murray and von Neumann proved that this is a total ordering: either $H _ { 0 } \leq H _ { 1 } ;$ or $H _ { 1 } \ \preceq \ H _ { 0 } ;$ or both, in which case $H _ { 0 }$ and $H _ { 1 }$ are isomorphic. For example, in a finite-dimensional type I situation, where H is a multiple of n copies of a single irreducible representation, each subrepresentation is the sum of $m \leqslant n$ copies of the irreducible representation, and the order structure of the (isomorphism classes of) subrepresentations is the same as the order structure of the integers $\{ 0 , 1 , \ldots , n \}$ .

Murray and von Neumann showed that the only order structures that can arise from factors are the following very simple ones:

$$
\text { Type   I }, \quad \{0, 1, 2, \dots , n \} \text {   or   } \{0, 1, 2, \dots , \infty \};
$$

$$
\text { Type   II }, \quad [ 0, 1 ] \text {   or   } [ 0, \infty ];
$$

$$
\text { Type   III }, \{0, \infty \}.
$$

The type of a factor is determined from the order structure of its projections according to this table.

In the case of factors of type II, the order structure is that of an interval of real numbers, not integers. Any subrepresentation of an isotypical representation of type II can be divided into yet smaller subrepresentations: we shall never reach an irreducible “atom.” Nevertheless, subrepresentations can still be compared in size by means of the “real-valued dimension” provided by Murray and von Neumann’s theorem.

A notable example of a factor of type II may be obtained as follows. Let G be a group and let $H = \ell ^ { 2 } ( G )$ be a Hilbert space having basis vectors [g] corresponding to the elements $g \in G$ . Then there is a natural representation of $G$ on H derived from the group multiplication law, called the regular representation: given an element $_ g$ of $G ,$ , the corresponding unitary map $U _ { g }$ is the linear operator that takes each basis vector $[ g ^ { \prime } ]$ in $\ell ^ { 2 } ( G )$ to the basis vector [gg	 ]. The commutant of this representation is a von Neumann algebra M. If G is a commutative group, then all the operators $U _ { g }$ are in the center of M; but if G is far enough from commutativity (for instance, if it is a free group), then M will have trivial center and will therefore be a factor. It can be shown that this factor is of type II. There is a simple explicit formula for the real-valued dimension of a subrepresentation corresponding to an orthogonal projection $P \in M .$ Represent P by an infinite matrix relative to the basis [g] of H. Because P commutes with the representation, it is easy to see that the diagonal elements of P are all the same, equal to some real number between 0 and 1. This real number is the dimension of the subrepresentation corresponding to P.

More recently, the Murray–von Neumann dimension theory has found unexpected applications in topology [I.3 §6.4]. Many important topological concepts, such as Betti numbers, are defined as the (integer-valued) dimensions of certain vector spaces. Using von Neumann algebras, one can define real-valued counterparts of these quantities that have useful additional properties. In this way, one can use von Neumann algebra theory to obtain topological conclusions. The von Neumann algebras used here are typically obtained by the construction of the previous paragraph from the fundamental group [IV.6 §2] of some compact space.

# 2.3 Modular Theory

Type III factors remained rather mysterious for a long time; indeed, Murray and von Neumann were at first unable to determine whether any such factors existed. They eventually managed to do so, but the fundamental breakthrough in the area came well after their pioneering work, when it was realized that each von Neumann algebra has a special family of symmetries, its so-called modular automorphism group.

To explain the origins of modular theory, let us consider once again the von Neumann algebra obtained from the regular representation of a group G. We defined the operators $U _ { g }$ on $\ell ^ { 2 } ( G )$ by multiplying on the left by elements of $G ;$ but we could equally well have considered a representation defined by multiplying on the right. This would have yielded a different von Neumann algebra.

So long as we deal only with discrete groups G this difference is unimportant, because the map $S : [ g ] \mapsto $ $[ g ^ { - 1 } ]$ ] is a unitary operator on H that interchanges the left and right regular representations. But for certain continuous groups the problem arises that the function $f ( g )$ may be square-integrable while $f ( g ^ { - 1 } )$ is not. In this situation there is no simple unitary isomorphism analogous to the one for discrete groups. To remedy this, one must introduce a correction factor called the modular function of $G .$

The project of modular theory is to show that something analogous to the modular function can be constructed for any von Neumann algebra. This object then serves as an invariant for all factors of type III, whether or not they are explicitly derived from groups.

Modular theory exploits a version of the GNS construction (section 1.4). Let M be a self-adjoint algebra of operators. A linear functional $\phi : M \to \mathbb { C }$ is called a state if it is positive in the sense that $\phi ( T ^ { * } T ) \geqslant 0 ,$ for every $T \in M$ (this terminology is derived from the connection described earlier between Hilbert space theory and quantum mechanics). For the purposes of modular theory we restrict attention to faithful states, those for which $\phi ( T ^ { * } T ) = 0$ implies $T = 0 .$ . If φ is a state, then the formula

$$
\langle T _ {1}, T _ {2} \rangle = \phi (T _ {1} ^ {*} T _ {2})
$$

defines an inner product on the vector space M. Applying the GNS procedure, we obtain a Hilbert space $H _ { M }$ . The first important fact about $H _ { M }$ is that every operator T in M determines an operator on $H _ { M }$ . Indeed, a vector $V \in H _ { M }$ is a limit $V = \operatorname* { l i m } _ { n \to \infty } V _ { n }$ of elements in $M ,$ and we can apply an operator $T \in M$ to the vector V using the formula

$$
T V = \lim _ {n \to \infty} T V _ {n},
$$

where on the right-hand side we use multiplication in the algebra M. Because of this observation, we can think of M as an algebra of operators on $H _ { M }$ , rather than as an algebra of operators on whatever Hilbert space we began with.

Next, the adjoint operation equips the Hilbert space $H _ { M }$ with a natural “antilinear” operator $S : H _ { M } \to H _ { M }$ by the formula1 $S ( V ) = V ^ { * }$ . Since $U _ { g } ^ { * } = U _ { g } .$ −1 for the regular representation, this is indeed analogous to the operator S we encountered in our discussion of continuous groups. The important theorem of Minoru Tomita and Masamichi Takesaki asserts that, as long as the original state $\phi$ satisfies a continuity condition, the complex powers $U _ { t } = ( S ^ { * } S ) ^ { \mathrm { i } t }$ have the property that $U _ { t } M U _ { - t } = M$ , for all t.

The transformations of M given by the formula $T \mapsto$ ${ U _ { t } T } U _ { - }$ t are called the modular automorphisms of M. Alain Connes proved that they depend only in a rather inessential way on the original faithful state $\phi .$ T o be precise, changing φ changes the modular automorphisms only by inner automorphisms, that is, transformations of the form $T \mapsto U T U ^ { - 1 }$ , where U is a unitary operator in M itself. The remarkable conclusion is that every von Neumann algebra M has a canonical oneparameter group of “outer automorphisms,” which is determined by M alone and not by the state φ that is used to define it.

The modular group of a type I or type II factor consists only of the identity transformation; however, the modular group of a type III factor is much more complex. For example, the set

$$
\{t \in \mathbb {R}: T \mapsto U _ {t} T U _ {- t} \text {   is   an   inner   automorphism } \}
$$

is a subgroup of R and an invariant of M that can be used to distinguish between uncountably many different type III factors.

# 2.4 Classification

A crowning achievement of von Neumann algebra theory is the classification of factors that are approximately finite dimensional. These are the factors that are in a certain sense limits of finite-dimensional algebras. Besides the range of the dimension function, which separates factors into types, the sole invariant is the module. This is a flow on a certain space that is assembled from the modular automorphism group.

A lot of attention is currently being given to the longstanding problem of distinguishing among the type II factors associated with the regular representations of groups. Of special interest is the case of free groups [IV.10 §2], around which has flourished the subject of free probability theory. Despite intensive effort, some fundamental questions remain open: at the time of writing it is unknown whether the factors associated with the free groups on two and on three generators are isomorphic.

Another important development has been subfactor theory, which attempts to classify the ways in which factors can be realized within other factors. A remarkable and surprising theorem of Vaughan Jones shows that, in the type II situation, where continuous values of dimensions are the norm, the dimensions of subfactors can in certain situations assume only a discrete range of values. The combinatorics associated with this result have also appeared in other apparently quite unrelated parts of mathematics, notably knot theory [III.44].

# 3 C∗-Algebras

Von Neumann algebra theory helps describe the structure of a single representation of a group on a Hilbert space. But in many situations it is of interest to gain an understanding of all possible unitary representations. To shed some light on this problem we turn to a related but different part of operator algebra theory.

Consider the collection B(H) of all bounded operators on a Hilbert space H. It has two very different structures: algebraic operations, such as addition, multiplication, and formation of adjoints; and analytic structures, such as the operator norm

$$
\| T \| = \sup \{\| T u \|: \| u \| \leqslant 1 \}.
$$

These structures are not independent of one another. Suppose, for instance, that T < 1 (an analytic hypothesis). Then the geometric series

$$
S = I + T + T ^ {2} + T ^ {3} + \dots
$$

converges in B(H), and its limit S satisfies

$$
S (I - T) = (I - T) S = I.
$$

It follows that I  T is invertible in (H) (an algebraic conclusion). One can easily deduce from this that the spectral radius r (T ) of any operator T (defined to be the greatest absolute value of any complex number in the spectrum of T ) is less than or equal to its norm.

The remarkable spectral radius formula goes much further in the same direction. It asserts that $r ( T ) ~ =$ $\scriptstyle \operatorname* { l i m } _ { n \to \infty } \| T ^ { n } \| ^ { 1 / n }$ . If T is normal $( T T ^ { * } \ = \ T ^ { * } T ) ,$ , and in particular if T is self-adjoint, then it may be shown that $\| T ^ { n } \| = \| T \| ^ { n }$ . As a result, the spectral radius of T is precisely equal to the norm of T . There is therefore a very close connection between the algebraic structure of (H), particularly algebraic structure related to the adjoint operation, and the analytic structure.

Not all the properties of (H) are relevant to this connection between algebra and analysis. $\mathrm { ~ A ~ } C ^ { * }$ -algebra A is an abstract structure that has enough properties for the argument of the previous two paragraphs to remain valid. A detailed definition would be out of place here, but it is worth mentioning that a crucial condition relating norm, multiplication, and -operation is

$$
\left\| a ^ {*} a \right\| = \left\| a \right\| ^ {2}, \quad a \in A,
$$

called the $C ^ { * }$ -identity for A. We also note that special classes of operators on Hilbert space (unitaries, orthogonal projections, and so on) all have their counterparts in a general C∗-algebra. For example, a unitary u ∈ A satisfies uu $\mathbf { \Psi } ^ { * } = u ^ { * } u = 1$ , and a projection p satisfies $p = p ^ { 2 } = p ^ { * }$ .

A simple example of a C∗-algebra is obtained by starting with a single operator $T \in { \mathcal { B } } ( H )$ . The collection of all operators $S \in { \mathcal { B } } ( H )$ that can be obtained as limits of polynomials in T and $T ^ { * }$ is ${ \mathrm { ~ \tt ~ { ~ \tt ~ } ~ } } _ { \mathrm { ~ \tt ~ { ~ C ~ } ~ } ^ { * } }$ -algebra said to be generated by T . The $C ^ { * }$ -algebra generated by T is commutative if and only if T is normal; this is one reason for the importance of normal operators.

# 3.1 Commutative $c ^ { * } { \bf \delta } { \bf - A }$ lgebras

If X is a compact [III.9] topological space [III.90], then the collection C(X) of continuous functions $f :$ X  C comes with natural algebraic operations (inherited from the usual ones on C) and a norm $\| f \| \bigstar =$ sup $\{ | f ( x ) | : x \in X \}$ . In fact, these operations make $C ( X )$ into a $C ^ { * }$ -algebra. The multiplication in C(X) is commutative, because the multiplication of complex numbers is commutative.

A basic result of Gelfand and Naimark asserts that every commutative $C ^ { * }$ -algebra is isomorphic to some C(X). Given a commutative $C ^ { * }$ -algebra A, one constructs X as the collection of all algebra homomorphisms $\xi : A \to \mathbb { C } ,$ , and the Gelfand transform then associates with $a \in A$ the function $\xi \mapsto \xi ( { a } )$ from X to C.

The Gelfand–Naimark theorem is a foundational result of operator theory. For example, a modern proof of the spectral theorem might proceed as follows. Let T be a self-adjoint or normal operator on a Hilbert space H, and let A be the commutative $C ^ { * }$ -algebra generated by T . By the Gelfand–Naimark theorem, A is isomorphic to C(X) for some space X, which may in fact be identified with the spectrum of T . If v is a unit vector in H, then the formula $S \mapsto \langle S \nu , \nu \rangle$ defines a state φ on A. The GNS space associated with this state is a Hilbert space of functions on $X ,$ and elements of $A \ = \ C ( X )$ act as multiplication operators. In particular, T acts as a multiplication operator. A small additional argument shows that T is unitarily equivalent to this multiplication operator, or at least to a direct sum of such operators (which is itself a multiplication operator on a larger space).

Continuous functions can be composed: if $f$ and $_ g$ are continuous functions (with the range of $_ g$ contained in the domain of $f ) ,$ then $f \circ g$ is also a continuous function. Since the Gelfand–Naimark theorem tells us that any self-adjoint element of a $C ^ { * }$ -algebra A sits inside an algebra isomorphic to the continuous functions on the spectrum of $^ { a , }$ we conclude that if $a \in A$ is selfadjoint, and if $\mathrm { ~ : ~ } f$ is a continuous function defined on the spectrum of $^ { a , }$ then an operator $f ( a )$ exists in A. This functional calculus is a key technical tool in $C ^ { * }$ -algebra theory. For example, suppose that $u \in A$ is unitary and $\| u - 1 \| < 2$ . Then the spectrum of u is a subset of the unit circle in C that does not contain 1. One can define a continuous branch of the complex logarithm function on such a subset, and it follows that there is an element a  log u of the algebra such that $a = - a ^ { * }$ and $u = \mathbf { e } ^ { a }$ . The path $t \mapsto \mathrm { e } ^ { t a } , 0 \leqslant t \leqslant 1$ , is then a continuous path of unitaries in A connecting u to the identity. Thus every unitary sufficiently close to the identity is connected to the identity by a unitary path.

# 3.2 Further Examples of $c ^ { * } .$ -Algebras

# 3.2.1 The Compact Operators

An operator on a Hilbert space has finite rank if its range is a finite-dimensional subspace. The operators of finite rank form an algebra, and its closure is a $C ^ { * } .$ algebra called the algebra of compact operators and denoted . One can also view  as a “limit” of matrix algebras

$$
M _ {1} (\mathbb {C}) \to M _ {2} (\mathbb {C}) \to M _ {3} (\mathbb {C}) \to \dots ,
$$

where each matrix algebra is included in the next by

$$
A \mapsto \left( \begin{array}{c c} A & 0 \\ 0 & 0 \end{array} \right).
$$

Many natural operators are compact, including the integral operators that arose in Fredholm’s theory. The identity operator on a Hilbert space is compact if and only if that Hilbert space is finite dimensional.

# 3.2.2 The CAR Algebra

The presentation of K as a limit of matrix algebras leads one to consider other “limi ${ [ s ^ { \prime } }$ of a similar sort. (We shall not attempt a formal definition of these limits here, but it is important to note that the limit of a sequence $A _ { 1 } \to A _ { 2 } \to A _ { 3 } \to \cdot \cdot \cdot$ depends on the homomorphisms $A _ { i } \ \to \ A _ { i + 1 }$ as well as on the algebras $A _ { i \cdot } )$ One particularly important example is obtained as the limit

$$
M _ {1} (\mathbb {C}) \to M _ {2} (\mathbb {C}) \to M _ {4} (\mathbb {C}) \to \dots ,
$$

where each matrix algebra is included in the next by

$$
A \mapsto \left( \begin{array}{c c} A & 0 \\ 0 & A \end{array} \right).
$$

This is called the CAR algebra, because it contains elements that represent the canonical anticommutation relations that arise in quantum theory. $C ^ { * }$ -algebras find several applications to quantum field theory and quantum statistical mechanics which extend von Neumann’s formulation of quantum theory in terms of Hilbert space.

# 3.2.3 Group C∗-Algebras

If G is a group and $g \mapsto U _ { g }$ is a unitary representation of G on a Hilbert space H, we can consider the smallest $C ^ { * }$ -algebra of operators on H containing all the $U _ { g } ;$ this is called the $C ^ { * }$ -algebra generated by the representation. An important example is the regular representation on the Hilbert space $\ell ^ { 2 } ( G )$ generated by $G ,$ which we defined in section $2 . 2 .$ . The $C ^ { * }$ -algebra that it generates is denoted $C _ { \mathrm { r } } ^ { * } \left( G \right)$ . The subscript $" \mathrm { r } "$ refers to the regular representation. Considering other representations leads to other, potentially different, group C∗-algebras.

Consider, for example, the case $G = \mathbb { Z }$ . Since this is a commutative group, its C∗-algebra is also commutative, and thus it is isomorphic to $C ( X )$ for a suitable $X ,$ by the Gelfand–Naimark theorem. In fact, X is the unit circle $S ^ { 1 }$ , and the isomorphism

$$
C (S ^ {1}) \cong C _ {\mathrm{r}} ^ {*} (\mathbb {Z})
$$

takes a function on the circle to its Fourier series.

States defined on group $C ^ { * }$ -algebras correspond to positive-definite functions defined on groups, and hence to unitary group representations. In this way new representations may be constructed and studied. For example, using states of group C∗-algebras it is possible to give to the set of irreducible representations of G the structure of a topological space.

# 3.2.4 The Irrational Rotation Algebra

The algebra $C ^ { * } ( \mathbb { Z } )$ is generated by a single unitary element U (corresponding to $1 \in \mathbb { Z } )$ . Moreover, it is the universal example of such a $C ^ { * }$ -algebra, which is to say that given any $C ^ { * }$ -algebra A and unitary u  A, there is one and only one homomorphism $C ^ { * } ( \mathbb { Z } ) \to A$ sending U to u. In fact, this is nothing other than the functional calculus homomorphism for the unitary u.

If instead we consider the universal example of a $C ^ { * } .$ algebra generated by two unitaries U, V subject to the relation

$$
U V = \mathrm{e} ^ {2 \pi \mathrm{i} \alpha} V U,
$$

where α is irrational, we obtain a noncommutative $C ^ { * }$ - algebra called the irrational rotation algebra $A _ { \alpha } .$ The irrational rotation algebras have been studied intensively from a number of points of view. Using K-theory (see below) it has been shown that $A _ { \alpha _ { 1 } }$ is isomorphic to $A _ { \alpha _ { 2 } }$ if and only if ${ \pm \alpha _ { 1 } \pm \alpha _ { 2 } }$ is an integer.

It can be shown that the irrational rotation algebra is simple, which implies that any pair of unitaries U, V satisfying the commutation relation above will generate a copy of $A _ { \alpha } .$ . (Note the contrast with the case of a single unitary: 1 is a unitary operator, but it does not generate a copy of $C ^ { * } ( \mathbb { Z } ) . )$ This allows us to give a concrete representation of $A _ { \alpha }$ on the Hilbert space $L ^ { 2 } ( S ^ { 1 } )$ , where U is the rotation through 2πα and V is multiplication by $z : S ^ { 1 } \to \mathbb { C }$ .

# 4 Fredholm Operators

A Fredholm operator between Hilbert spaces is defined to be a bounded operator T for which the kernel and cokernel are finite dimensional. This means that the homogeneous equation $T u \ = \ 0$ admits only finitely many linearly independent solutions, while the inhomogeneous equation $T u = \nu$ admits a solution if v satisfies a finite number of linear conditions. The terminology arises from Fredholm’s original work on integral equations; he showed that if K is an integral operator, then $I + K$ is a Fredholm operator.

For the operators that Fredholm considered, the dimensions of the kernel and cokernel must be equal, but in general this need not be so. The unilateral shift operator S, which maps the infinite “row vec-$\mathrm { t o r " } \left( a _ { 1 } , a _ { 2 } , a _ { 3 } , \ldots \right)$ to $( 0 , a _ { 1 } , a _ { 2 } , \dots )$ , is an example. The equation $S u \ : = \ : 0$ has only the zero solution, but the equation $S u \ = \ \nu$ has a solution only if the first coordinate of the vector v is zero.

The index of a Fredholm operator is defined to be the integer difference

$$
\operatorname{index} (T) = \dim (\ker (T)) - \dim (\operatorname{coker} (T)).
$$

For example, every invertible operator is a Fredholm operator of index 0, whereas the unilateral shift is a Fredholm operator of index 1.

# 4.1 Atkinson’s Theorem

Consider the two systems of linear equations

$$
\left\{ \begin{array}{l} 2. 1 x + y = 0 \\ 4 x + 2 y = 0 \end{array} \right\} \quad \text { and } \quad \left\{ \begin{array}{l} 2 x + y = 0 \\ 4 x + 2 y = 0 \end{array} \right\}.
$$

Although the coefficients of these equations are very close, the dimensions of their kernels are quite different: the left-hand system has only the zero solution, whereas the right-hand system has the nontrivial solutions $( t , - 2 t )$ . Thus the dimension of the kernel is an unstable invariant of the system of equations. A similar remark applies to the dimension of the cokernel. By contrast, the index is stable, despite its definition as the difference of two unstable quantities.

An important theorem of Frederick Atkinson gives precise expression to these stability properties. Atkinson’s theorem asserts that an operator T is Fredholm if and only if it is invertible modulo compact operators. This implies that any operator that is sufficiently close to a Fredholm operator is itself a Fredholm operator with the same index, and that if T is a Fredholm operator and K is a compact operator, then $T + K$ is a Fredholm operator with the same index as T . Notice that, since integral operators are compact operators, this contains Fredholm’s original theorem as a special case.

# 4.2 The Toeplitz Index Theorem

topology [I.3 §6.4] studies those properties of mathematical systems that remain the same when the system is (continuously) perturbed. Atkinson’s theorem tells us that the Fredholm index is a topological quantity. In many contexts it is possible to obtain a formula for the index of a Fredholm operator in terms of other, apparently quite different, topological quantities. Formulas of this sort often indicate deep connections between analysis and topology and often have powerful applications.

The simplest example involves the Toeplitz operators. A Toeplitz operator has a matrix with the special form

$$
T = \left( \begin{array}{c c c c c} b _ {0} & b _ {1} & b _ {2} & b _ {3} & \dots \\ b _ {- 1} & b _ {0} & b _ {1} & b _ {2} & \dots \\ b _ {- 2} & b _ {- 1} & b _ {0} & b _ {1} & \dots \\ b _ {- 3} & b _ {- 2} & b _ {- 1} & b _ {0} & \dots \\ \vdots & \vdots & \vdots & \vdots & \ddots \end{array} \right).
$$

In other words, as you go down each diagonal of the matrix, the entries remain constant. The sequence of coefficients $\{ b _ { n } \} _ { n = - \infty } ^ { \infty }$ defines a function $f ( z ) ~ =$ $\scriptstyle \sum _ { n = - \infty } ^ { \infty } b _ { n } z ^ { - n }$ =−∞ on the unit circle in the complex plane, called the symbol of the Toeplitz operator. It can be shown that a Toeplitz operator whose symbol is a continuous function which is never zero is Fredholm. What is its index?

The answer is given by thinking about the symbol as a mapping from the unit circle to the nonzero complex numbers: in other words, as a closed path in the nonzero complex plane. The fundamental topological invariant of such a path is its winding number: the number of times it “goes around” the origin in the counterclockwise direction. It can be proved that the index of a Toeplitz operator with nonzero symbol f is minus the winding number of f . For example, if f is the function $f ( z ) ~ = ~ z$ (with winding number 1), then the associated Toeplitz operator is the unilateral shift S that we encountered earlier (with index 1). The Toeplitz index theorem is a very special case of the atiyah–singer index theorem [V.2], which gives a topological formula for the indices of various Fredholm operators that arise in geometry.

# 4.3 Essentially Normal Operators

Atkinson’s theorem suggests that compact perturbations of an operator are in some sense “small.” This leads to the study of properties of an operator that are preserved by compact perturbation. For instance, the essential spectrum of an operator T is the set of complex numbers λ for which T  λI fails to be Fredholm (that is, invertible modulo compact operators). Two operators $T _ { 1 }$ and $T _ { 2 }$ are essentially equivalent if there is a unitary operator U such that $U T _ { 1 } U ^ { * }$ and $T _ { 2 }$ differ by a compact operator. A beautiful theorem originally due to weyl [VI.80] asserts that two self-adjoint or normal operators are essentially equivalent if and only if they have the same essential spectrum.

One might argue that the restriction to normal operators in this theorem is inappropriate. Since we are concerned with properties that are preserved by compact perturbation, would it not be more appropriate to consider essentially normal operators—that is, operators T for which $T ^ { * } T - T T ^ { * }$ is compact? This apparently modest variation leads to an unexpected result. The unilateral shift S is an example of an essentially normal operator. Its essential spectrum is the unit circle, as is the essential spectrum of its adjoint; however, S and $S ^ { * }$ cannot be essentially equivalent, because S has index 1 and $S ^ { * }$ has index 1. Thus some new ingredient, beyond the essential spectrum, is needed to classify essentially normal operators. In fact, it follows easily from Atkinson’s theorem that if essentially normal operators $T _ { 1 }$ and $T _ { 2 }$ are to be essentially equivalent, then not only must they have the same essential spectrum but also, for every λ not in the essential spectrum, the Fredholm index of $T _ { 1 }$ λI must be equal to the Fredholm index of $T _ { 2 } - \lambda I .$ . The converse of this statement was proved by Larry Brown, Ron Douglas, and Peter Fillmore in the 1970s, using entirely novel techniques that led to a new era of interaction between $C ^ { * }$ -algebra theory and topology.

# 4.4 K-Theory

A remarkable feature of the Brown–Douglas–Fillmore work was the appearance within it of tools from algebraic topology [IV.6], notably K-theory. Remember that, according to the Gelfand–Naimark theorem, the study of (suitable) topological spaces and the study of commutative $C ^ { * }$ -algebras are one and the same; all the techniques of topology can be transferred, via the Gelfand–Naimark isomorphism, to commutative $C ^ { * } .$ algebras. Having made this observation, it is natural to ask which of these techniques can be extended further, to provide information about all $C ^ { * }$ -algebras, commutative or not. The first and best example is K-theory.

In its most basic form, K-theory associates with each C∗-algebra A an Abelian group $K ( A )$ , and with each homomorphism of $C ^ { * }$ -algebras a corresponding homomorphism of Abelian groups. The building blocks for $K ( A )$ can be thought of as generalized Fredholm operators associated with $A ;$ the generalization is that these operators act on “Hilbert spaces” in which the complex scalars are replaced by elements of the $C ^ { * }$ -algebra $A .$ The group $K ( A )$ itself is defined to be the collection of connected components of the space of all such generalized Fredholm operators. Thus if $A = \mathbb { C } ,$ for instance (so that we are dealing with classical Fredholm operators), then $K ( A ) = \mathbb { Z } .$ . This follows from the fact that two Fredholm operators are connected by a path of Fredholm operators if and only if they have the same index.

One of the great strengths of K-theory is that one can construct K-theory classes from a variety of different ingredients. For example, every projection $p \in A$ defines a class in $K ( A )$ which can be thought of as a “dimension” for the range of $p .$ This connects K-theory to the classification of factors (section 2.2), and has become an important tool in the effort to classify various families of $C ^ { * }$ -algebras, such as the irrational rotation algebras. (It was at one time thought that the irrational rotation algebras might not contain any nontrivial projections at all: the construction of such projections by Marc Rieffel was an important step in the development of $C ^ { * }$ -algebra K-theory.) Another beautiful example is George Elliott’s classification theorem for locally finite-dimensional $C ^ { * }$ -algebras like the CAR algebra; they are completely determined by K-theoretic invariants.

The problem of computing the K-theory groups of noncommutative $C ^ { * }$ -algebras, particularly group $C ^ { * } \cdot$ algebras, has turned out to have important connections with topology. In fact, some key advances in topology have come from $C ^ { * }$ -algebra theory in this way, thereby allowing operator algebraists to repay some of the debt they owe to the topologists for K-theory. The principal organizing problem in this area is the Baum– Connes conjecture, which proposes a description of the K-theory of group $C ^ { * }$ -algebras in terms of invariants familiar in algebraic topology. Most of the progress on the conjecture to date is the result of work of Gennadi Kasparov, who dramatically broadened the original discoveries of Brown, Douglas, and Fillmore to cover not just single essentially normal operators but also noncommuting systems of operators, that is, $C ^ { * }$ -algebras. Kasparov’s work is now a central component of operator algebra theory.

# 5 Noncommutative Geometry

descartes’s [VI.11] invention of coordinates showed that one can do geometry by thinking about coordinate functions rather than directly thinking about points in space and their interrelationships: these coordinate functions are the familiar x, y, and z. The Gelfand– Naimark theorem can be viewed as one expression of this idea of passing from the “point picture” of a space $X$ to the “field picture” of the algebra $C ( X )$ of functions on it. The success of K-theory in operator algebras invites us to ponder whether the field picture might be more powerful than the point picture, since K-theory can be applied to noncommutative $C ^ { * }$ -algebras which may not have any “points” (homomorphisms to C) at all.

One of the most exciting research frontiers in operator algebra theory is reached along a path which develops these thoughts. The noncommutative geometry program of Connes takes seriously the idea that a general $C ^ { * }$ -algebra should be thought of as an algebra of functions on a “noncommutative space,” and goes on to develop “noncommutative” versions of many ideas from geometry and topology, as well as completely new constructions that have no commutative counterpart. Noncommutative geometry begins with the creative reformulation of ideas from ordinary geometry in ways that involve only operators and functions, but not points.

Consider, for instance, the circle $S ^ { 1 }$ . The algebra $C ( S ^ { 1 } )$ ) reflects all the topological properties of $S ^ { 1 }$ , but to incorporate its metric (distance-related) properties as well we look not just at $C ( S ^ { 1 } )$ but at the pair consisting of the algebra $C ( S ^ { 1 } )$ and the operator $D = \mathrm { i d } / \mathrm { d } \theta$ on the Hilbert space $H = L ^ { 2 } ( S ^ { 1 } )$ . Notice that if f is a function on the circle (considered as a multiplication operator on $H ) ,$ , then the commutator Df $f D$ is also a multiplication operator, this time by i df /dθ. It follows that ordinary measurements of angular distance between points on the circle can be recovered from $C ( S ^ { 1 } )$ and D by the formula

$$
d (p, q) = \max \left\{\left| f (p) - f (q) \right|: \| D f - f D \| \leqslant 1 \right\}.
$$

Connes argues that operator $| D | ^ { - 1 }$ plays the role of the “unit of arc-length ${ \mathrm { d } } s "$ in this and many other, more complicated situations.2

Another feature of the examples Connes considers, also of central importance in noncommutative geometry, is the fact that the operator $| D | ^ { - k }$ is a traceclass operator (see section 1.5) when k is large enough. In the case of the circle, k needs to be bigger than 1. Computations with traces connect noncommutative geometry to cohomology theory [IV.6 §4]. We now have two kinds of “noncommutative algebraic topology,” namely K-theory and a new variant of homology called cyclic cohomology; the connection between the two is provided by a very general index theorem.

There are several procedures that produce noncommutative C∗-algebras (to which Connes’s methods can be applied) from classical geometric data. The irrational rotation algebras $A _ { \theta }$ are examples; the classical picture to which they apply is the quotient space [I.3 §3.3] of the circle by the group of rotations through multiples of θ. Classical methods of geometry and topology are unable to handle this quotient space, but the noncommutative approach via $A _ { \theta }$ is much more successful.

An exciting but speculative possibility is that the basic laws of physics should be addressed from the perspective of noncommutative geometry. The transition to noncommutative $C ^ { * }$ -algebras can be viewed as analogous to the transition from classical to quantum mechanics. However, Connes has argued that noncommutative $C ^ { * }$ -algebras play a role in describing the physical world even before the transition is made to quantum physics.

# Further Reading

Connes, A. 1995. Noncommutative Geometry. Boston, MA: Academic Press.   
Davidson, K. 1996. C∗-Algebras by Example. Providence, RI: American Mathematical Society.   
Fillmore, P. 1996. A User’s Guide to Operator Algebras. Canadian Mathematical Society Series of Monographs and Advanced Texts. New York: John Wiley.   
Halmos, P. R. 1963. What does the spectral theorem say? American Mathematical Monthly 70:241–47.   
2. The operator D is not quite invertible since it vanishes on constant functions. A small modification must therefore be made before considering inverse operators. The operator D is by definition the positive square root of $D ^ { 2 } .$ .

# IV.16 Mirror Symmetry

# Eric Zaslow

# 1 What Is Mirror Symmetry?

Mirror symmetry is a phenomenon found in theoretical physics that has had profound mathematical applications. It burst onto the mathematical scene after Candelas, de la Ossa, Green, and Parkes exploited the physical phenomenon to make precise predictions about certain sequences of numbers describing geometric spaces. The sequence predicted by those authors began 2875, 609 250, 317 206 375, . . . , and was far beyond the scope of calculation at the time. The phenomenon of mirror symmetry is that some physical theories have equivalent, “mirror” theories that lead to the same predictions. If some prediction requires a hard calculation but is easy to perform in the mirror theory, then you can get the answer for free! These physical theories do not have to be realistic models of physics. For instance, beginning students of physics often study point particles on frictionless planes. Although they are unrealistic, such toy models can bring the physical concepts into focus and their analysis can give rise to very interesting mathematics.

# 1.1 Exploiting Equivalences

Children at school in the 1950s used log tables to exploit the equivalence of multiplication of positive numbers with addition of real numbers. Given the problem of multiplying two large numbers a and $^ { b , }$ they would use a table to look up the logarithms log(a) and log(b) (to a certain number of significant figures), then add them by hand. They would then use the same table to find which number had a logarithm equal to log(a)  log(b). The answer is ab.

College students sometimes exploit the equivalence defined by fourier transforms [III.27] to solve differential equations. Basically, the Fourier transform is a rule that maps one function $f ( x )$ to a new function ${ \hat { f } } ( { \boldsymbol { p } } )$ . What is nice is that the transform of the derivative $f ^ { \prime } ( x )$ relates in a very simple way to ${ \hat { f } } ( { \boldsymbol { p } } )$ : it is $\mathrm { i } p \hat { f } ( p )$ , where i is the imaginary number $\sqrt { - 1 }$ . If you want to solve a differential equation such as $f ^ { \prime } ( x ) + 2 f ( x ) = h ( x )$ , where h(x) is a given function and you are trying to find f , you can map the equation to its Fourier transform equation $\mathrm { i } p { \hat { f } } ( p ) + 2 { \hat { f } } ( p ) =$ $\hat { h } ( p )$ . This is much easier: it is an algebraic equation rather than a differential equation, and has the solution $\hat { f } ( p ) = \hat { h } ( p ) / ( 2 + \mathrm { i } p )$ . The solution f (x) is then the function which has $\hat { h } ( p ) / ( 2 + \mathrm { i } p )$ as its Fourier transform.

Mirror symmetry is like a fancy Fourier transform, mapping much more information than is contained in a single function. Every aspect of a physical theory is involved.

This article will (eventually) focus on the mathematics of mirror symmetry, but it is crucial to understand its physical origins. We therefore begin with a brief guide to physics. (For a further discussion of mathematical physics, see vertex operator algebras [IV.17 §2].) This is in no way an adequate treatment—a separate Companion to Physics would be needed—but we hope to give enough of the flavor of the subject to help the reader with the later sections. (A reader familiar with physical theories may wish to skip the next section and refer back as needed.)

# 2 Theories of Physics

# 2.1 Formulations of Mechanics and Action Principles

# 2.1.1 Newtonian Physics

Newton’s second law states that a particle moving through space accelerates1 in proportion to the force it experiences: $F = m { \ddot { x } }$ . The force is itself the (negative) gradient of a gravitational potential V (x), so this equation can be written mx¨ $+ ~ nabla V ( x ) = 0$ . Stationary particles sit at minima of the potential: examples are a ball in equilibrium at the end of a spring, or a pea at the bottom of a bowl. In stable situations, there is a restoring force proportional to some displacement distance. This means that in some appropriate coordinate, $F \sim - x ,$ , so $V ( x ) = k x ^ { 2 } / 2$ , for some k. The solutions are oscillatory, with angular frequency $\omega = \sqrt { k / m }$ . This model is called the simple harmonic oscillator.

# 2.1.2 The Least Action Principle

Every major theory can also be formulated by means of an idea known as the least action principle. Let us see how it works for the equations of Newtonian mechanics. Consider an arbitrary path of a particle x(t) and form the quantity

$$
S (x) = \int \left[ \frac {1}{2} m \dot {x} ^ {2} - V (x) \right] d t.
$$

Here and below, the notation x may represent more than one coordinate. If x is used as a point in spacetime, it will include the time coordinate, if that is not otherwise noted. Likewise, we omit component notation on most vectors. The notation should be clear from the context. The quantity S(x), which is known as the action, equals the kinetic energy minus the potential energy. One then considers which paths minimize this action. That is, we ask which paths x(t) have the property that, when they are perturbed by a small amount $\delta x ( t )$ , the action is unchanged, to leading order. (So in fact we require only that the action is unchanged to first order, and not that it is actually minimized. Solutions of saddle-point type are allowed.) The answer turns out to be precisely those paths that satisfy $m { \ddot { x } } + \nabla V ( x ) = 0 . ^ { 2 }$

For example, consider the simple harmonic oscillator in two dimensions. We can model x as a complex number and set $V ( x ) = k | x | ^ { 2 }$ . The action is then $\int { \frac { 1 } { 2 } } [ m | \dot { x } | ^ { 2 } - k | x | ^ { 2 } ]$ . Note that a phase rotation $x \to \mathrm { e } ^ { \mathrm { i } \theta } x$ leaves the action invariant, and is therefore a symmetry of the equations of motion.

Lesson. Physical solutions extremize the action.

The principle of least action applies to many other physical situations, as we shall see below. First, though, we describe another formulation of mechanics.

# 2.1.3 The Hamiltonian Formulation of Mechanics

hamilton’s [VI.37] formulation of the equations of motion also deserves mention. It leads to first-order equations. Let S be the action and define L by $S = \int L { \mathrm { d } } t ,$ , and consider the (typical) case where L is a function of coordinates x and their time derivatives x˙. Then set $p = \mathrm { d } L / \mathrm { d } \dot { x }$ , a function that can depend both on x and on x˙. (In the example $L = { \textstyle { \frac { 1 } { 2 } } } m { \dot { x } } ^ { 2 } - V ( x )$ that we have already considered, we find that $p = m { \dot { x } } ,$ or $\dot { x } = p / m . )$ Now let us consider the function $H = p { \dot { x } } - L$ , which is called the hamiltonian [III.35], and change variables from (x, x)˙ to (x, p) so as to remove all mention of x˙. In the example, H works out to be

$$
\frac {p ^ {2}}{m} - \left(\frac {p ^ {2}}{2 m} - V (x)\right) = \frac {p ^ {2}}{2 m} + V (x),
$$

which is the total energy. For the simple harmonic oscillator, $H = p ^ { 2 } / 2 m + k x ^ { 2 } / 2$ .

The equations $\dot { x } = \partial H / \partial p$ and $\dot { p } = -  { \partial } H /  { \partial } x$ are the equations of motion in the Hamiltonian formulation; they can be shown to be equivalent to those obtained from the action principle. In the example, ${ \dot { x } } = p /$ m and $\dot { p } = - \nabla V$ . Using the first equation to replace p by mx˙ in the second, we recover the equation $m { \ddot { x } } + \nabla V ( x ) =$ 0. More generally, one can consider the time derivative of some quantity $f ( x , p )$ constructed from p and x and prove—using the chain rule and the equations of motion—that

$$
\dot {f} = \frac {\partial f}{\partial x} \frac {\partial H}{\partial p} - \frac {\partial f}{\partial p} \frac {\partial H}{\partial x} = \{H, f \}.
$$

The term in the middle is called the Poisson bracket of H and $f ,$ denoted $\{ H , f \}$ .

Lesson. The Hamiltonian controls time dependence through the Poisson bracket.

Notice that when we plug the coordinates x and p themselves into the bracket, we derive the identity

$$
\{x, p \} = - 1. \tag {1}
$$

It is also possible to begin with the Hamiltonian viewpoint. One considers a space endowed with a bracket operation on functions, such that there are coordinate functions (not uniquely determined) obeying $\{ x , p \} =$ −1. The mechanical model is defined by a function $H ( x , p )$ , which determines the dynamics.

# 2.1.4 Symmetry

A brief remark on symmetry is in order. noether [VI.76] proved that in the action formulation of mechanics, a symmetry of the action results in a conserved quantity. The prototypical example is translational or rotational symmetry, where the potential of a particle is invariant under some direction of translation or rotation: the corresponding conserved quantity is then momentum or angular momentum. In the example above, $V ( x ) = k | x | ^ { 2 } / 2$ is independent of $\theta ,$ the phase of x. The equation of motion determined by varying θ is $\mathsf { d } ( m | x | ^ { 2 } \dot { \theta } ) / \mathsf { d } t \ = \ 0$ , so in this case it is the angular momentum $m | x | ^ { 2 } { \dot { \theta } }$ that is conserved. In the Hamiltonian formulation, since a conserved quantity $f ( x , p )$ does not change with time, it must have zero Poisson bracket with the Hamiltonian: $\{ H , f \} = 0 \quad$ . In particular, the Hamiltonian itself is conserved.

# 2.1.5 Action Functions for Other Theories

Returning now to action principles, we shall see how different physical theories are described through different actions. In electricity and magnetism, maxwell’s equations [IV.13 §1.1] can be formulated in the form $\delta S = 0 ,$ where now the action S takes the form of an integral over space and time of the electric (E) and magnetic (B) fields. In the case where there are no sources, the action is written

$$
S = \frac {1}{8 \pi e ^ {2}} \int [ E ^ {2} - B ^ {2} ]   \mathrm{d} x   \mathrm{d} t, \tag {2}
$$

where e is the electric charge of an electron. There is one important difference from the previous example, which is that the variations of the action must be taken with respect to the fundamental fields, and E and B are not fundamental as they are derived from the electromagnetic potential $A = ( \phi , A )$ by the equations $E \ = \ \nabla \phi - { \dot { A } } , \ B \ = \ \nabla \times { \boldsymbol { A } } .$ If you rewrite S in terms of A, vary A by δA, and set $\delta S = 0$ , then you recover Maxwell’s equations from the least action principle.

It is clear that the electromagnetic action merely changes sign under the replacement ${ \cal E }  { \cal B } , { \cal B }  - { \cal E } ,$ and therefore any solution $\delta { \cal S } ~ = ~ 0$ remains a solution under the transformation. This is an example of an equivalence of a classical theory of physics. In fact, this symmetry extends to the case where there are sources (such as electrons) if we also interchange electric and magnetic sources. (No magnetic sources have been observed in the universe, but a theory with such objects still makes sense.)

Lesson. Physical equivalences act on fields and their sources.

Electricity and magnetism is a “field theory,” which means that the degrees of freedom involve functions that depend on position in space. Contrast this with Newtonian mechanics, where the spatial degrees of freedom are just the coordinates of the particle(s). However, there is not much conceptual distance between the two, as can be seen in the following toy model.

We will consider the simplest example: a scalar field, $\phi .$ That is, φ is just a function that takes numerical values. Now imagine that space has just one dimension, not three, and further that that dimension is a circle, which we can describe with an angular coordinate, θ. At any fixed point in time we can use fourier series [III.27] to write the scalar field as $\phi ( \theta ) \quad = \quad$ $\textstyle \sum _ { n } c _ { n } \exp ( \mathrm { i } n \theta )$ , where the $c _ { n }$ are the Fourier coefficients, and if we want the values of φ to be real numbers then we must insist that $c _ { - n } = c _ { n } ^ { * }$ . We can then think of $\phi ( \theta )$ not as a function but as an infinitedimensional vector $( c _ { 0 } , c _ { 1 } , \ldots )$ . The spatial dependence of $\phi$ is completely determined by the coefficients $c _ { n } .$ . If we now wish to consider time dependence, then all we have to do is use time-dependent components $( c _ { 0 } ( t ) , c _ { 1 } ( t ) , \ldots )$ , which looks a lot like an infinite set of quantum-mechanical particles $c _ { n } .$ Thus, the function φ has the Fourier expansion $\begin{array} { r } { \phi ( \theta , t ) = \sum _ { n } c _ { n } ( t ) } \end{array}$ exp(inθ).

The simplest action for a scalar field $\phi$ that allows wave-like solutions of the equations of motion serves as a natural analogue of equation (2):

$$
S = \int \frac {1}{2 \pi} [ (\dot {\phi}) ^ {2} - (\phi^ {\prime}) ^ {2} ] d \theta d t, \tag {3}
$$

where $\phi ^ { \prime } = \partial \phi / \partial \theta$ . When we plug the Fourier expansion into the action and perform the θ integration, we get

$$
S = \int \sum_ {n} [ | \dot {c} _ {n} | ^ {2} - n ^ {2} | c _ {n} | ^ {2} ] d t. \tag {4}
$$

Note that the term in brackets is just the action for a particle $c _ { n }$ in a quadratic potential, as in section 2.1.2. We simply have an infinite number of harmonic oscillators (with the exception of the $c _ { 0 }$ degree of freedom, which corresponds to a free particle in no potential).

Lesson. Field theory is like point particle theory with an infinite number of particles. The particles correspond to the degrees of freedom of the field. When the action is just quadratic in the derivatives, the particles have an interpretation as simple harmonic oscillators.

Even general relativity [IV.13] fits into this framework as a field theory. For a space-time $M ,$ the field is the riemannian metric [I.3 §6.10] on space-time. The metric is what determines the lengths of paths between points—so a stretching of space-time, for example, is represented by a rescaled metric. The action is then constructed as the integral of the Riemannian curvature scalar over space-time: $\begin{array} { r } { S = \int _ { M } \mathcal { R } . } \end{array}$ 3

# 2.2 Quantum Theory

Mirror symmetry is an equivalence of quantum theories, so we must develop an understanding of what a quantum theory is and what an equivalence looks like. There are two formulations of quantum mechanics: the operator formulation and Feynman’s path-integral formulation.

Both formulations are probabilistic, meaning that you cannot predict exactly what will be observed in a single measurement, but you can make precise predictions about what will be observed after multiple, repeated measurements in the same environment. For instance, your experimental apparatus may involve a beam of electrons hitting a screen and making a mark. The beam will contain millions of electrons, so the pattern of marks on the screen can be predicted with great accuracy. However, we cannot say what will happen to a single, given electron—all we can do is assign probabilities to the outcomes of various measurements. These probabilities are encoded in the so-called “wave function” Ψ of the particle.

# 2.2.1 Hamiltonian Formulation

In the operator formulation of quantum mechanics, the positions and momenta of classical mechanics (and any quantity formed from them) are converted into operators [III.50] acting on a hilbert space [III.37] according to the following rule: replace the Poisson bracket $\{ \cdot , \cdot \} b y i / \hbar [ \cdot , \cdot ]$ , where $[ A , B ] = A B - B A$ is the commutator bracket and  is Planck’s constant. Thus, for example, we get from equation (1) the relation $[ x , p ] = \mathrm { i } \hbar .$ . The state of a particle (or system) is now defined not as a set of values of x and p but as a vector Ψ in the Hilbert space. Once again, time evolution is determined by the Hamiltonian, H, but now H is an operator. The basic dynamical equation is

$$
H \Psi = \mathrm{i} \hbar \frac {\mathrm{d}}{\mathrm{d} t} \Psi . \tag {5}
$$

This is called the Schrödinger equation.

Lesson. To quantize a classical theory, replace ordinary degrees of freedom by operators on a vector space; replace Poisson brackets by commutator brackets.

In the case where we have a particle on the real line R, the Hilbert space is the space of square-integrable functions $L ^ { 2 } ( \mathbb { R } )$ , so we write Ψ as $\Psi ( x )$ . The commutation relation is obeyed if we think of x as the operator that sends the function Ψ (x) to the function xΨ (x). Now the relation $[ x , p ] = \mathrm { i } \hbar$ means that we should represent p as the operator i(d/dx). The values of the classical quantity associated with an operator correspond to the eigenvalues [I.3 §4.3] of that operator, so for example a state with momentum p has the form $\cal { Y } ~ \sim ~ \exp ( \mathrm { i } p x / \hbar )$ . Unfortunately, this is not squareintegrable on the real line, but it would become so if we identified x and $x + 2 \pi R .$ , for some number (radius) $R \ > \ 0$ . Topologically, this compactifies [III.9] R to a circle, but note that Ψ will be single-valued only if $p = n \hbar / R _ { ☉ }$ , where n is an integer. Thus, momentum is “quantized” in units of $\hbar / R . ^ { 4 }$ The integer label of the $c _ { n }$ of equation (4) can therefore also be thought of as a momentum.

In the above example, R is the degree of freedom of the classical coordinate x. In other examples, there is a copy of $L ^ { 2 } ( \mathbb { R } )$ for each real degree of freedom, whether or not it represents a geometric location.

Another novelty is that position and momentum do not commute as operators in quantum mechanics, meaning they cannot be simultaneously diagonalized: you cannot specify the position and momentum simultaneously. This is a form of Heisenberg’s uncertainty principle (see operator algebras [IV.15 §1.3]).

# 2.2.2 Symmetry

As the rules of quantization would suggest, a symmetry of a quantum theory is an operator A such that $[ H , A ] = 0$ . That ${ \mathrm { i } } s ,$ A commutes with the Hamiltonian, and therefore respects the dynamics.

# 2.2.3 Example: The Simple Harmonic Oscillator

We now discuss an example that will be useful later on for understanding quantum field theory and mirror symmetry: the simple harmonic oscillator in quantum mechanics. Suppose that the constants are chosen so that the Hamiltonian is given by $H = x ^ { 2 } + p ^ { 2 }$ . If one defines $a = \left( x + \mathrm { i } p \right) / \sqrt { 2 }$ and $a ^ { \dagger } = ( x - \mathrm { i } p ) / \sqrt { 2 }$ , then one can show that $a ^ { \dagger }$ raises the energy of a state by one unit5 and a lowers the energy by one unit. Invoking the physical argument that there is a ground state $\psi _ { 0 }$ of lowest energy, this state must obey $a \psi _ { 0 } = 0 .$ . One then finds that all states can be written in terms of the basis vectors $\Psi _ { n } = ( a ^ { \dagger } ) ^ { n } \Psi _ { 0 }$ with energy $n + { \frac { 1 } { 2 } }$ . Note that $\psi _ { 0 }$ has energy $\frac { 1 } { 2 } . ^ { 6 }$ The basis $\{ \psi _ { n } \}$ is called the occupation number basis, since the interpretation is that $\Psi _ { n }$ has n energy “quanta” above the ground state.

# 2.2.4 Path-Integral Formulation

Feynman’s path integral formulation of quantum mechanics builds on the idea of the least action principle. In this formulation, the probability of an experiment is calculated through an average over all paths of particles, and not just the ones which extremize the action. Each path x(t) is weighted by the factor exp(iS(x)/), where S(x) is the action of the path x(t) and  is Planck’s constant, which is very small compared with macroscopic action scales. This average can be an imaginary number, but the probability of the process is the square of its absolute value.

Note that $\exp ( \mathrm { i } S / \hbar ) = \cos ( S / \hbar ) + \mathrm { i } \sin ( S / \hbar )$ , so if S changes appreciably when we vary x(t), then the real and imaginary parts will oscillate rapidly, since  is small. Then, when we integrate over paths x(t), the positive and negative oscillations will roughly cancel. As a result, the main contributions to the weighted sum over paths will come from those paths for which S does not vary when the path does: the classical paths! However, if the variations are sufficiently small compared with $\hbar ,$ then nonclassical paths can contribute appreciably. One typically separates the degrees of freedom into the classical trajectory piece and the quantum fluctuations near it. Then one can organize the path integral in a perturbation theory around the parameter .

We have not yet discussed the integrand of the path integral, and will not go into the details of this. The main point is that the theory makes a prediction about the likelihood of measuring a physical process. Each process determines a possible integrand. For example, from our discussion above we learn that the integrand for measuring the likelihood of a quantum-mechanical particle going from the point $x _ { 0 }$ at time $t _ { 0 }$ to the point $x _ { 1 }$ at time $t _ { 1 }$ gives nonzero weight—determined by the exponentiated action—to all paths that go from $x _ { 0 }$ to $x _ { 1 }$ as t goes from $t _ { 0 }$ to $t _ { 1 }$ and zero weight to all other paths.

It is illustrative to consider a toy model of a path integral on a “space-time” that consists of just a single point. Then the possible “paths” of a scalar field, say, are simply the values that the field can take at the point, so they are real numbers. The action is then an ordinary function S(x) on R. For the purposes of this example, let us consider the case where $\mathrm { i } S / \hbar = - \frac { 1 } { 2 } x ^ { 2 } + \lambda x ^ { 3 }$ . The possible integrands are (sums of) powers of x, so the basic path integrals to perform are $\int x ^ { k } \exp ( - { \frac { 1 } { 2 } } x ^ { 2 } +$ $\lambda x ^ { 3 } )$ dx, which we denote by $\langle x ^ { k } \rangle$ . The value at $\lambda = 0$ is easily calculated.7 For small λ we expand $\mathrm { e } ^ { \lambda x ^ { 3 } }$ as $1 +$ $\lambda x ^ { 3 } + \lambda ^ { 2 } x ^ { 6 } / 2 + \cdot \cdot \cdot$ , and evaluate each term by the same methods as for $\lambda = 0$ . This is how we construct a welldefined perturbation theory, even when the integral is not calculable.

As we see from this example, path integrals are easiest when the action is only quadratic in the variables, just as we found in the operator formulation of quantum mechanics. The mathematical reason for this is that Gaussian integrals (exponentials of squares) can be done explicitly, while integrals involving exponentials of cubics or higher are difficult or impossible. For quadratic actions, the path integral can be evaluated exactly, but when cubic or higher terms appear, the perturbation series is necessary.

# 2.2.5 Quantum Field Theory

The generalization to field theories follows our earlier pattern. We think of quantum field theories, then, as being like quantum mechanics with infinite numbers of particles. In fact, the quantum field theories in which the fields Φ and their derivatives do not have more than quadratic terms in the action are easily understood in this way—we had a preview of this in equation (4). The Fourier components correspond to particles indexed by their momenta. Each one looks like a simple harmonic oscillator at some frequency, which will depend on the Fourier coefficient. The quantum Hilbert space is then a (tensor) product of lots of different “occupation number Hilbert spaces,” one for each Fourier component of each field. Since the occupation number basis is also an energy eigenbasis, these states have a simple time evolution under the Hamiltonian H. That is, if $H = E$ on some state $\Psi ( t = 0 )$ , then that state evolves like

$$
\Psi (t) = \exp (\mathrm{i} E t / \hbar) \Psi (0).
$$

However, if the action includes terms that are cubic or higher, then things get interesting: particles can decay! This can be seen, for example, from the scalar field of equation (3) if we include a term $\phi ^ { 3 }$ in the action, and therefore also the Hamiltonian. If we write this using Fourier components, we get terms involving three oscillators, such as $a _ { 3 } ^ { \dag } a _ { 4 } ^ { \dag } a _ { 7 }$ . To see this, recall that after we quantize the real field φ, the Fourier components

$c _ { n }$ act as harmonic oscillators, and we have written $a _ { n }$ for the associated creation and annihilation operators. Since the Hamiltonian governs time evolution according to equation (5), this means that over time one particle (the 7 mode) can decay into two others (the 3 and the 4). Such decay processes occur in real life, and it is a great triumph of quantum field theory that it can predict such events with astounding accuracy.

In fact, because the space of paths of fields is infinite dimensional, the path integral in quantum field theory is not usually defined in a mathematically rigorous way. However, the perturbation series for producing predictions can be defined just as for quantum mechanics, and this is how physicists make their predictions in practice. This perturbation series is organized in terms of Feynman diagrams (which are discussed in vertex operator algebras [IV.17]). These diagrams, and the rules for computing them, completely solve the perturbation problem.

As in the example of quantum mechanics, different integrands of the path integral correspond to different predictions. If Φ is some function of the fields of some quantum field theory, we write Φ for the path integral with Φ as an integrand (as we did for $\langle x ^ { k } \rangle$ in the previous section). We call such a term a “correlation function.” If $\phi = \phi _ { 1 } ( x _ { 1 } ) \cdot \cdot \cdot \phi _ { n } ( x _ { n } )$ , the answer will depend on the action of the theory, the fields $\phi _ { i } ,$ , and the space-time points $x _ { i }$

One might wonder if a symmetry of a classical theory always remains a symmetry of the same theory after quantization. The answer is sometimes no. Such a case is known as an “anomaly.” Roughly speaking, this is because the measure of integration of the path integral is not preserved under the symmetry, but this is a somewhat heuristic explanation because the path integral has no rigorous definition in general.

Returning to our cubic example, if the interaction term $\phi ^ { 3 }$ has a coefficient λ, so that it is $\lambda \phi ^ { 3 }$ , then we organize the perturbation series as a power series in λ. In terms of paths, probabilities of decay processes can be evaluated by considering paths that split into two— like the letter Y—with each leg carrying the label of the appropriate particle.

# 2.2.6 String Theory

Feynman’s perturbation theory has an important generalization in string theory. String theory considers particles not as points but as loops. Instead of paths of particles through space-time, we get paths of loops,

$$
\begin{array}{l} \int \exp (- \frac {1}{2} x ^ {2} + J x) d x = \int \exp (- \frac {1}{2} (x + J) ^ {2}) \exp (\frac {1}{2} J ^ {2}) d x \\ = \sqrt {2 \pi} \exp (\frac {1}{2} J ^ {2}). \\ \end{array}
$$

Now if we differentiate this answer with respect to J, and set $J = 0 ,$ , we get x . Taking k derivatives gives xk , and the theory is solved.

which look like two-dimensional surfaces. String theory amplitudes are computed by summing over all surfaces. These sums are organized in a perturbation series in powers of the so-called string coupling constant, $\lambda _ { g }$ . The power of $\lambda _ { g }$ in the perturbation series depends on the number of holes in the surface.

The surfaces are called worldsheets. At each point of the worldsheet, its location in space-time is determined by coordinates $X ^ { i }$ . These coordinates themselves depend on the location on the worldsheet. In effect, we get an auxiliary theory: a field theory of coordinates on the two-dimensional surface! In string theory, even this two-dimensional field theory must be considered as a quantum field theory. The fields of the two-dimensional theory are maps from the surface to actual space-time. However, from the point of view of the worldsheet, the worldsheet itself is a twodimensional space-time and the maps are fields on this space-time with values in some other (target) space.

Mirror symmetry was discovered as a result of the study of these quantum field theories on two-dimensional surfaces. Subsequently, the same phenomenon was discovered in the case where the strings were not closed loops but filaments with endpoints. Both cases play an important role below.

# 3 Equivalence in Physics

Mirror symmetry is a particular type of equivalence of quantum field theories. As we have seen, quantum field theories are rules for producing probabilities of physical processes. In the path-integral formulation, probabilities are computed from correlation functions of fields. According to Feynman, these correlation functions can be thought of as being averages over all paths of fields. Each path is weighted by exp(iS/), where S is the action of the path and  is Planck’s constant. Let us denote the correlation function of some integrand Φ in theory A as $\langle \phi \rangle _ { \mathrm { ~ , ~ } }$ A. Recall that Φ can depend on various fields $\phi _ { i }$ and points of space-time $x _ { i } ,$ , and the correlation function will depend on all these and the action of theory A.

Equivalence, then, is a map from all possible fields $\phi _ { i }$ in a theory A to corresponding fields $\tilde { \phi } _ { i }$ in a theory B such that

$$
\langle \Phi \rangle_ {A} = \langle \tilde {\Phi} \rangle_ {B}.
$$

(For the moment, we deliberately neglect to notate the dependence on the points $x _ { i \cdot } )$ One special correlation function is 1 , which we call the partition function and denote by Z. As the field 1 always gets mapped to 1, we derive the corollary that the partition functions must be equal: $Z _ { \mathrm { A } } = Z _ { \mathrm { B } }$ .

Of course, this all has a description in the operator formulation of the quantum theory. Each state Ψ and each operator a in one theory must get mapped to a corresponding state $\tilde { \psi }$ and operator a˜ in the mirror theory, in such a way that corresponding operators map corresponding states to states which themselves correspond. Here one sees the sharp analogy with the slide rule and the operations of multiplication and addition of numbers.

Each theory is typically described through some mathematical model, so an equivalence implies a host of mathematical identities between quantities constructed from corresponding models.

The particular case of mirror symmetry refers to an equivalence of quantum field theories on a two-dimensional surface. The most typical example of mirror symmetry is the physical theory whose fields are maps ϕ from a two-dimensional riemann surface [III.79] Σ to some target space, M. Such a theory is called a sigma model. As we saw above, in string theory M plays the role of actual space-time, but for our purposes we can even consider the case where M is the real line $\mathbb { R } ,$ s o that ϕ is an ordinary function. This case has already been studied in section 2.1.5. The action is given in equation (4). We can then write the partition function as

$$
Z = \langle 1 \rangle = \int [ \mathcal {D} \varphi ] \mathrm{e} ^ {\mathrm{i} S (\varphi) / \hbar},
$$

where [Dϕ] represents the measure of integration over all paths.8

One approach to evaluating the partition function Z is through a process known as Wick rotation. One first Euclideanizes the time coordinate by writing $\tau \ = \ \mathrm { i } t$ (this is the Wick rotation), which leads to an imaginary Euclidean action $\mathrm { i } S _ { \mathrm { E } } .$ . One then tries to evaluate the path integral in this framework, hoping that the answer will be holomorphic [I.3 §5.6]. If it is, then one can use analytic continuation to work out the answer for ordinary time. The advantage is that the Euclidean exponential weighting becomes exp( SE/), so the minima of $S _ { \mathrm { E } }$ receive the greatest weighting and the integral might converge. The nonconstant minima of the Euclidean action are called instantons. After Euclideanizing equation (4), the action becomes the “energy” SE of the

map ϕ:

$$
S _ {\mathrm{E}} = \int_ {\Sigma} | \nabla \varphi | ^ {2}.
$$

The energy of a map has a conformal symmetry, meaning that it is independent of local scale transformations on the Riemann surface, that is, transformations that can be locally approximated by a combination of rotations and dilations. Invariance under rescaling by a positive number λ can easily be seen: each of the two derivatives in $| \nabla \varphi | ^ { 2 }$ decreases by a factor of λ, while the area element increases by $\lambda ^ { 2 } .$ . Rotational invariance is clear from the form of $| \nabla \varphi | ^ { 2 }$ . The combination of the two, along with the fact that this argument did not depend on the derivatives of the scaling parameter λ, leads to the statement of local scale invariance.

The conformal symmetry of the action is an example of a classical symmetry of the action that is not necessarily maintained in the quantum theory. However, the quantum theory has no anomaly—meaning that the symmetry is preserved—if M is chosen to be a complex, calabi–yau manifold [III.6].

The Calabi–Yau condition can be thought of as a complex notion of orientation. Recall that for an oriented manifold one can continuously choose, on each patch, a basis for the tangent space such that, when we move from patch to patch, the determinant of the changeof-basis matrix is equal to one. The same is true on a Calabi–Yau manifold, but now we consider complex bases for the complex tangent spaces.

When the target manifold is a Calabi–Yau manifold, the instantons are complex analytic maps from the twodimensional surface. Instantons are not “close” to the constant paths; their effects are therefore not accessible by perturbative methods such as Feynman diagrams. They are therefore “nonperturbative” phenomena. An example from quantum mechanics would be a particle in a double-well potential such as $( x ^ { 2 } - 1 ) ^ { 2 }$ . The zero-energy minima are the two constant (stationary) paths at x 1. An instanton path could go from $x = - 1 \mathrm { t o } x = + 1$ , or vice versa. Such trajectories occur and are known as “quantum tunneling.”

Lesson. Inaccessible by perturbation theory, instantonic effects are notoriously challenging to calculate.

# 3.1 Mirror Pairs

In the setting above, we considered maps from a twodimensional surface Σ to a target (Calabi–Yau) space. Let us denote this quantum field theory by Q(M), which is shorthand for the collection of all fields and all possible correlation functions created from them. In this setup, we say that the Calabi–Yau manifolds M and W are “mirror pairs” if Q(M) is equivalent to $Q ( W )$ . Through the magic of mirror symmetry, hard problems in Q(M) involving instantons can be answered in Q(W ) by considering only the much simpler constant paths.

# 4 Mathematical Distillation

A physical theory contains a tremendous amount of information. For example, correlation functions can involve any number of fields, each evaluated at different points on the two-dimensional surface. This is typically too unwieldy a situation to approach mathematically. Instead, equipped with a symmetry of the theory called “supersymmetry,” a mathematical distillation can be performed. The distillation procedure is called topological twisting, and the resulting “topological field theory” has correlation functions that are independent of the positions of points. Because of this independence, the correlation functions are certain characteristic numbers associated with the underlying geometric setup. In fact, there are two types of twisting, typically called A and B, which capture different aspects of the manifold in question.

# 4.1 Complex and Symplectic Geometry

# 4.1.1 Complex Geometry

To get a feel for the geometric aspect captured by topological twisting, recall that we can construct the circle $S ^ { 1 }$ from the real line R by identifying the points θ and $\theta + 2 \pi$ , and therefore also $\theta + 2 \pi n$ , where n is any integer. What we have done is identified points related by a lattice of integer translations. We could choose the lattice to consist of multiples of some other real number $r ,$ but since any two such lattices differ only by an overall scaling of R, we would effectively get the same space. In the complex plane C, we can do the same thing with a two-dimensional lattice of translations generated by two complex numbers $\lambda _ { 1 }$ and $\lambda _ { 2 } ,$ as long as the quotient $\lambda _ { 2 } / \lambda _ { 1 }$ is not real. This space is called a torus and has the same topology as any two-dimensional surface with one hole. It has more structure, however, because it can be covered by regions described by a complex coordinate—with different regions related by complex analytic maps. The pairs $( \lambda _ { 1 } , \lambda _ { 2 } )$ and $( \lambda _ { 1 } , \lambda _ { 2 } ~ + ~ \lambda _ { 1 } )$ generate the same lattice of translations, as do the pairs $( \lambda _ { 1 } , \lambda _ { 2 } )$ and $( \lambda _ { 2 } , - \lambda _ { 1 } )$ . In fact, lattices related by a complex rescaling of C are equivalent, so a better parametrization of the lattice is the ratio $\tau = \lambda _ { 2 } / \lambda _ { 1 }$ .

By redefining the direction of one of the λs, we can assume that the imaginary part of τ is positive, so τ takes values in the upper half of the complex plane. By the reasoning above, we note that τ and $\tau + 1$ , as well $\mathrm { a s } - 1 / \tau$ , all come from the same lattice. The number τ can also be thought of in the following way. The torus has two distinct loops, one generated by a straight path from $z \ \mathrm { t o } \ z + \lambda _ { 1 }$ , and one generated by a straight path from $z \ \mathrm { t o } \ z + \lambda _ { 2 }$ . Then $\lambda _ { 1 }$ and $\lambda _ { 2 }$ are both the result of the line integral of the complex differential dz over the loop. In fact, the loop did not even need to be straight to lead to this conclusion. The values of such integrals over subspaces without boundaries (the loops, here) are more generally called periods.

Although any two tori are topologically equivalent, one can show that there is no complex analytic map between two complex tori described by genuinely different values of τ. The parameter τ therefore determines the complex geometry of the space. Roughly speaking, we think of this parameter as describing the shape of the torus. (See moduli spaces [IV.8 §2.1] for a further discussion of this.)

The topological B-model depends only on the complex geometry of the target space M. That is, the theory depends, continuously, only on the parameter τ.

# 4.1.2 Symplectic Geometry

Another aspect of geometry is the size of the torus, which is described simply by an area element. Let us recall that, topologically, all tori look like $\mathbb { R } ^ { 2 }$ with points identified by the lattice of integer horizontal and vertical translations (but not necessarily in a way that would respect any complex geometry). The points of the torus can be thought of as the unit square with opposite sides glued together. An area element in $\mathbb { R } ^ { 2 }$ looks like ρ dx dy, which then determines the area $\rho$ of the unit square. These notions of two-dimensional area generalize to two-dimensional subspaces in higher-dimensional spaces. The study of such structures is called symplectic geometry [III.88], and so we call $\rho$ the symplectic parameter.

The topological A-model depends only on the symplectic geometry of the target space M. That is, the theory depends, continuously, only on the parameter $\rho .$

# 4.2 Cohomological Theories

As you might imagine, the passage from an ordinary theory to a topological theory involves identifying many aspects of the physical theory that were previously distinct, such as different point values of a single field. Mathematically, a well-established method of producing topological aspects of a structure—and one that involves making identifications—is through a cohomology theory [IV.6 §4]. Cohomology theories follow the pattern of having an operator δ obeying the equation $\delta \circ \delta = 0$ . We think of this equation as the statement image $( \delta ) \subset \ker ( \delta )$ . The cohomology group H(δ) is formed as the quotient $H ( \delta ) ~ = ~ \ker ( \delta ) / \operatorname* { i m a g e } ( \delta )$ , which means that we identify any two vectors u and v satisfying $\delta u \ = \ \delta \nu \ = \ 0$ , so long as the difference $u - \upsilon$ can be written as δw for some w. Then $H ( \delta )$ is just the space of all such vectors, up to identifications.

The topological twisting of physical theories is similar. The operator δ is a physical operator acting on a Hilbert space of states. The presence of supersymmetry in our theories ensures that δ exists and squares to zero. The vector states of the topological theory are just the elements of $H ( \delta )$ , i.e., states in the original theory Ψ obeying $\delta \Psi = 0 ,$ up to identification. In many cases, these states can be identified with ground states.

It is crucial that supersymmetry is a symmetry that contains the complex translations of points on the twodimensional surface. This means that the value of a field operator $\phi ( z )$ at one point is identified with its value $\phi ( z ^ { \prime } )$ at another. In other words, the physics of the topological theory is independent of the positions of the operators! In the path-integral formulation, this means that the correlation functions are independent of the positions of the fields inserted into the integrand. What can they depend on, then? They depend on the particular field or combination of fields inserted, and they depend on the geometric parameter (such as $\rho$ or τ) of the space M.

# 4.2.1 The A-Model and the B-Model

Given a Calabi–Yau space, one can actually construct two operators, $\delta _ { \mathrm { { A } } }$ and $\delta _ { \mathrm { B } } .$ , each of which squares to zero. There are therefore two distinct corresponding topological twistings and two distinct topological theories that can be constructed from a Calabi–Yau space.

If M and W are mirror Calabi–Yau pairs, you might wonder if the topological models constructed from them will still be equivalent theories. The answer is a most interesting form of yes: the resulting A-model of one Calabi–Yau manifold M is equivalent to the B-model of the mirror W , and vice versa! The complex and symplectic aspects of the theories get interchanged under mirror symmetry! In particular, a hard symplectic question of M might get mapped to an easy computation involving the complex geometry of W .

We emphasize here that the two manifolds may be completely topologically distinct. For example, the Euler characteristic of one is the negative of the other.

# 5 Basic Example: T-Duality

Although the circle is not complex, it provides a very illustrative entry into mirror symmetry that can be studied quite easily. We will find an equivalence between two theories constructed from circles. The equivalence will be very nontrivial, however, as states of very different kinds will be shown to correspond.

Consider the case where the two-dimensional surface is a cylinder, with spatial dimension a unit circle, and one dimension of time, and let us look at the sigma model (these were introduced in section 3). Suppose also that the target space is a circle of radius R, which we denote by $S _ { R } ^ { 1 }$ . We think of $S _ { R } ^ { 1 }$ as the real line, with two points identified if they differ by a multiple of 2πR. Maps from one circle to another can be classified by their winding number, an integer that tells you how many (net) times the image of a point goes around the second circle when the point goes once around the first. The map $\theta \mapsto m R \theta$ from the circle to $S _ { R } ^ { 1 }$ has winding number m. This allows us to write the field $\varphi ( \theta )$ as a winding piece, mRθ, plus an honest Fourier series (no winding): $\textstyle \varphi ( \theta ) = m R \theta + x + \sum _ { n \neq 0 } c _ { n }$ exp(inθ). Here we have singled out the constant mode $x \ : = \ : c _ { 0 }$ of the Fourier series. We have expanded just the θ dependence in a series, so every continuous parameter (x and the $c _ { n } )$ should be thought of as a function of time, as well.

The energy, or Hamiltonian, of such a map is computed as in section 2.1.3:

$$
H = (m R) ^ {2} + \dot {x} ^ {2} + \sum_ {n} | \dot {c} _ {n} | ^ {2} + n ^ {2} | c _ {n} | ^ {2}.
$$

Comparing this with the harmonic oscillator Hamiltonian of section 2.1.3, we can see that each degree of freedom $c _ { n } ( t )$ plays the role of a (complex) quantummechanical particle in a simple harmonic oscillator potential. There is an occupation-mode basis for describing the quantum mechanics of each mode.9 The full Hilbert space of the quantum theory is the (tensor) product of each of these, plus parts involving the constant mode and winding number, which we now discuss. (Remember, each degree of freedom of the classical theory becomes a particle in the quantum field theory.)

The constant mode x has energy $\dot { x } ^ { 2 }$ , and therefore has no associated potential (it can be anywhere on the circle). This mode represents a free quantum-mechanical particle on the circle. Recall that the momentum of the x particle is represented by the operator $- \mathrm { i } ( \mathrm { d } / \mathrm { d } x )$ . This operator has eigenfunctions $\mathrm { e } ^ { \mathrm { i } p x }$ . The requirement that these eigenfunctions are invariant under the translation $x \to x + 2 \pi R$ means that the eigenvalues of momentum are “quantized,” and have the form $p = n / R$ .

In contrast to momentum, the integer winding number (m) is really a classical label for the possible maps from a circle to a circle. Although integral, it is clearly on a different footing from the integer n of momentum. Still, it is also an important label on the Hilbert space. For each m, we have a space of m-winding configurations which gets quantized to become the mth sector of the Hilbert space. Roughly, this sector ${ \mathcal { H } } _ { m }$ comprises the functions of all the degrees of freedom of all the m-winding maps. We can consider the winding number as an operator by simply declaring that the states with winding number m have eigenvalue mR.

Ignoring the oscillator modes for the moment, the state of momentum $n / R$ with winding m has energy $( n / R ) ^ { 2 } + ( m R ) ^ { 2 }$ . In particular, the energy is unchanged if we make the simultaneous switches $( m , n )  ( n , m )$ and $R \  \ 1 / R .$ Since the oscillator modes $a _ { n }$ have energies that are independent of $R ,$ and since the modes are noninteracting particles, this symmetry can be extended to a full equivalence of the theories with targets $S _ { R } ^ { 1 }$ and $S _ { 1 / R } ^ { 1 } ,$ , with momentum in one theory corresponding to winding number in the other.

In this example, the target space $S ^ { 1 }$ is neither complex nor symplectic. As a result, we cannot construct the topological A- and B-models. Nevertheless, we have demonstrated the stronger statement that the two sigma models with target space $S _ { R } ^ { 1 }$ and $S _ { 1 / R } ^ { 1 }$ are equivalent. The theories are mirror pairs. In the special case of circles, mirror symmetry is referred to as T-duality. In fact, the entire phenomenon of mirror symmetry—even for noncircles—can be deduced from T-duality.

# 5.1 Tori

If we take the product of two circles $S _ { R _ { 1 } } ^ { 1 } \times S _ { R _ { 2 } } ^ { 1 }$ , we get a torus. We can think of the torus as a circle family of circles, since for each point in $S _ { R _ { 2 } } ^ { 1 }$ we have a circle $S _ { R _ { 1 } } ^ { 1 }$ . As we have seen in section 4.1.1, this space is complex— specifically, it is the complex plane C quotiented by a lattice of translations. A particularly simple lattice is the one generated by the translations $z  z + R _ { 1 }$ and $z  z + \mathrm { i } R _ { 2 }$ . As discussed in section 4.1.1 above, the lattice is determined by the complex number $\tau = \mathrm { i } R _ { 2 } / R _ { 1 }$ , equal to the ratio of integrals (“periods”) of the complex form dz over the two nontrivial loops of the torus.

The symplectic data is captured by the area element. Recall that we can choose coordinates x and $_ y$ such that the identifications look like unit translations in each direction. Then the (normalized) area element of the torus with radii $R _ { 1 }$ and $R _ { 2 }$ is $R _ { 1 } R _ { 2 }$ dx ${ \mathrm { d } } y$ , which integrates to $R _ { 1 } R _ { 2 }$ on the unit square. Let us define the symplectic parameter $\rho = \mathrm { i } R _ { 1 } R _ { 2 }$ . We now perform T-duality for the first circle $R _ { 1 } ~  ~ 1 / R _ { 1 }$ . We see that under this substitution, the complex and symplectic parameters get interchanged:10

$$
\tau \longleftrightarrow \rho .
$$

Lessons. Mirror symmetry interchanges complex and symplectic parameters. Mirror symmetry is T-duality.

# 5.2 The General Case

The torus is the only compact one-dimensional Calabi– Yau space and is therefore the simplest one, but the discussion above is part of a more general picture. The Calabi–Yau condition ensures a unique complex volume element, or orientation (dz, above), whose “periods” determine, and in turn vary with, the complex parameters. Though the A- and B-models both turn out to be rather simple in the case of the torus, what is important in general is that the B-model is completely determined by how the periods of the complex volume element (which were $\lambda _ { 1 }$ and λ2 in section 4.1.1) change with the parameters of the theory (of which there was just one in section 4.1.1, namely τ). Again, the relation $\tau = \lambda _ { 2 } / \lambda _ { 1 }$ is quite simple for the torus, but more complicated in general. In any case, this data gives all the information of the B-model. The reason for all of this is that the instantons of the B-model turn out to be just the constant maps. Each point of the target space determines a constant map, and as a result the B-model is reduced to (classical) complex geometry of the target space. This is determined by the periods.

This state of affairs is to be compared with the Amodel. The A-model depends on the symplectic parameters $\rho ,$ i.e., the areas of two-dimensional surfaces inside the target space. In contrast to the B-model, however, the dependence on ρ is very complicated, in general. The reason for this is that the instantons of the A-model are area-minimizing surfaces inside the target space, and their enumeration is a notoriously challenging problem. (The problem is not terribly challenging for the torus, however.) Mathematically, the A-model instantons are described by the theory of Gromov–Witten invariants, the subject to which we now turn.

# 6 Mirror Symmetry and Gromov–Witten Theory

As we mentioned above, the B-model on W is explained entirely by the classical complex geometry of W . The only relevant maps for B-model computations are the constant ones, so the space of such maps is equal to W itself, and correlators reduce to (classical) integrals over W . In fact, one of the integrands to be integrated is the complex volume element. Let us call the parameter for all possible complex volume elements τ. B-model correlation functions are then determined by τ-dependent integrals over W . In particular, the partition function $Z _ { \mathrm { B } } ^ { ( W ) }$ of the B-model on W depends on τ, so we write it as Z(W )B (τ). $Z _ { \mathrm { B } } ^ { ( W ) } ( \tau )$

The main point about topological twisting is that local variations of the fields are all identified, as they are related by the operator δ. In particular, varying the point on the worldsheet is a trivial operation in the topological theory. It turns out that, for the B-model on W , only the constant maps contributed, but for the A-model the situation is a bit more subtle. To give a feel for the geometry, consider again the winding of a map from a circle to a circle. Maps with different windings can never be deformed continuously into one another. The winding number is a measure of how the first circle “wraps” (or winds) around the target, according to the map. Because it is a discrete parameter it cannot change under continuous variations. Likewise, when M is a higher-dimensional space, the two-dimensional surface Σ can “wrap” around two-dimensional subspaces of M by different amounts. The parameters for wrapping are again discrete. A map ϕ can wrap Σ around the basic surfaces $C _ { i }$ in M by different integer amounts, $k _ { i }$ . We say that $k \ = \ k _ { i }$ labels the “class” of the map $\varphi .$ (More precisely, ϕ(Σ) is a closed 2-cycle when Σ is compact, and k labels its homology class.) Different classes k contribute through different (Euclidean) actions $S _ { k } ( \rho )$ , which depend on the areas $\rho$ and the class k but not on the continuous details of the map $\varphi _ { k }$ . The partition function can have contributions from all classes. Different classes may contribute differently not only through the exponential weighting, but also in accordance with how many minimal surfaces they contain. (A good example of a minimal surface in threedimensional space is a soap film. If you fix the boundary with a wire, the soap film will seek to find the minimumarea surface with that boundary.) In our examples, the space M is actually complex; the minimal surfaces we speak of in Gromov–Witten theory are complex analytic maps from Σ. That is, if you have a complex coordinate for Σ, then the complex coordinates for the surfaces M can be written as complex analytic functions of Σ.

The difference between the A-model and the B-model comes from the fact that the topological model is constructed from an operator $\delta ,$ which was guaranteed to exist by the presence of supersymmetry in our theories. For the different models, the relevant supersymmetry operators $\delta _ { \mathrm { { A } } }$ and $\delta _ { \mathrm { B } }$ are simply different. As we saw above, the maps relevant to the A-model are the instantons, or complex analytic maps from Σ to M. Roughly, then, A-model correlation functions on M, and in particular the partition function $Z _ { \mathrm { A } } ^ { ( { \cal M } ) }$ , are sums over classes k of surfaces in M and sums over instantons in each class, each one weighted by its instanton action $\exp ( - S _ { k } ( \rho ) )$ ). We have explicitly written the dependence on the parameter for the symplectic structure $\rho .$ . For Calabi–Yau manifolds, such maps should be discrete, and it is a conjecture, true in all known cases, that they are finite in number if we fix the class, k. All this data is packaged in a function of $\rho ,$ , and based on what we have argued, the partition function must take the general form

$$
Z _ {\mathrm{A}} ^ {(M)} (\rho) = \sum_ {k} n _ {k} \exp (- S _ {k} (\rho)).
$$

The coefficients $n _ { k }$ are called Gromov–Witten invariants.11 $a n t s . ^ { 1 1 }$

Putting things together, if $( M , A )$ is mirror to $( W , B )$ , and if we can identify for each complex parameter τ for W a corresponding symplectic parameter $\rho ( \tau )$ for M, then we have

$$
Z _ {\mathrm{A}} ^ {(M)} (\rho) = Z _ {\mathrm{A}} ^ {(M)} (\rho (\tau)) = Z _ {\mathrm{B}} ^ {(W)} (\tau). \tag {6}
$$

The first equality means we should rewrite $\rho$ in terms of τ, and the second says that the answer should be given by the corresponding B-model on W . Therefore, all of the information about complex analytic surfaces in M, which is encapsulated in the coefficients $n _ { k } ,$ is completely determined by the classical geometry of W !

This remarkable predictive power—the computation of an infinite number of difficult Gromov–Witten invariants through equations such as (6)—is what led to such intense interest in mirror symmetry at its inception.

# 7 Orbifolds and Nongeometric Phases

# 7.1 Nongeometric Theories

Mirror symmetry is about an equivalence of quantum field theories, and not every such field theory has the geometric content of a target space as in the sigma model. The structure involved in mirror symmetry—or at least its topological version—begins with a quantum theory with a supersymmetry algebra that allows for the passage to a topological theory. That is, there is a Hilbert space of states, a Hamiltonian operator, and a particular algebra of symmetries, i.e., operators that commute with the Hamiltonian. There are no dictates as to how one constructs such a setup, and the sigma model of maps to a target space is only one such way. Other methods abound. The geometric case is merely the one most suited for mathematicization (and exposition), which is why we have focused on the theory with a target space.

As an intermediate case—possibly geometric, possibly not—we will discuss the so-called orbifold theories.

# 7.2 Orbifolds

When space-time is a cylinder $S ^ { 1 } \times \mathbb { R } ,$ with a circle $S ^ { 1 }$ as its spatial dimension, there is a fascinating construction in quantum field theory known as an orbifold theory. This is defined as follows. Suppose there is a finite group G of symmetries (such as a reflection symmetry). That ${ \mathrm { i } } s ,$ each group element acts as an operator on the Hilbert space, so if $g \in G$ then it sends a state $\boldsymbol { \Psi }$ to a state $g \Psi .$ . Then one defines a new theory by identifying states related by the symmetry. To construct the theory, let us first consider the ground state $\psi _ { 0 }$ of the original theory. This is assumed to be invariant under the group: that is, $g \psi _ { 0 } = \Psi _ { 0 }$ for all group elements, $g . { } ^ { 1 2 }$ One then constructs the space ${ \mathcal { H } } _ { 0 }$ of all invariant states. This is known as the untwisted sector, and $\psi _ { 0 }$ is the ground state of the untwisted sector. In the case where G is commutative, a twisted sector is then constructed for every group element $g \in G . ^ { 1 3 }$ To construct the twisted sector, first think of the spatial dimension $S ^ { 1 }$ as being an interval [0, 1] with endpoints 0 and 1 identified. Recall that the Hilbert space of states is constructed from (functions of) all the degrees of freedom of the possible configurations of fields. The twisted sector ${ \mathcal { H } } _ { g }$ corresponds to additional field configurations Φ that are related at the two ends by the action of $g \colon \ s 0 \ \phi ( 1 ) \ = \ g \phi ( 0 )$ . Such field configurations represent configurations on the circle $S ^ { 1 }$ since left and right ends are related by the group, and therefore get identified. These additional configurations are thus part of the orbifold theory. One constructs a sector ${ \mathcal { H } } _ { g }$ of the Hilbert space by taking all such states $\Psi _ { g }$ that also obey the invariance condition $h \Psi _ { g } \ = \ \Psi _ { g }$ for all group elements h.

Orbifolds may be geometric, as they are in the case of the sigma model to a manifold X on which a discrete group G acts. For example, rotations act on the plane, and we can consider the four-element group generated by a right-angle rotation. The quotient of the plane by these rotations looks like a cone. As another example, the finite groups of symmetries of the platonic solids (tetrahedron, cube, etc.) act on the twodimensional sphere by rotations. When we take $X = S ^ { 2 }$ and G a platonic group, we get an interesting orbifold. In fact, if we simply take the space of orbits of the group G, it is topologically just a sphere again, but not a smooth one—it has cone points. These cone points would be troublesome in a quantum field theory, but the “stringy” orbifold is perfectly “smooth.”

The orbifold theory itself carries a symmetry. For example, if G is the commutative group with two elements, then there is an untwisted sector and a unique twisted sector. There is a symmetry corresponding to multiplication by 1 in the untwisted sector and by 1 in the twisted sector. This symmetry is not geometric. Orbifold theories with symmetries can often themselves be orbifolded in such a way as to recover the original theory. In fact, the theory and its orbifold are also often mirror pairs! Greene and Plesser used such a construction to create the first examples of mirror pairs. Furthermore, they used ways of ascribing geometric interpretations to some nongeometrically constructed theories so as to identify mirror Calabi–Yau spaces. To be precise, they took the space of all nonzero complex 5-vectors $X ~ = ~ ( X _ { 1 } , X _ { 2 } , X _ { 3 } , X _ { 4 } , X _ { 5 } )$ satisfying the equation

$$
X _ {1} ^ {5} + X _ {2} ^ {5} + X _ {3} ^ {5} + X _ {4} ^ {5} + X _ {5} ^ {5} + \tau X _ {1} X _ {2} X _ {3} X _ {4} X _ {5} = 0,
$$

identifying X with λX for any nonzero complex number λ. (If X is a solution, then so is λX.) The equation actually defines a family of complex spaces, since $\tau \in \mathbb { C }$ is a parameter. The orbifold theory is defined from the finite group of phase transformations

$$
\begin{array}{l} (X _ {1}, X _ {2}, X _ {3}, X _ {4}, X _ {5}) \\ \mapsto (\omega^ {n _ {1}} X _ {1}, \omega^ {n _ {2}} X _ {2}, \omega^ {n _ {3}} X _ {3}, \omega^ {n _ {4}} X _ {4}, \omega^ {n _ {5}} X _ {5}), \\ \end{array}
$$

where $\omega = \mathrm { e } ^ { 2 \pi \mathrm { i } / 5 }$ and $\textstyle \sum _ { i = 1 } ^ { 5 } n _ { i }$ is a multiple of 5. This space and its orbifold are actually the mirror pair about which Candelas et al. made their famous predictions.

# 8 Boundaries and Categories

The entire story of mirror symmetry becomes much richer when we allow the strings to have endpoints. Strings with ends are called “open strings,” while “closed strings” refers to loops. Mathematically, allowing ends corresponds to adding boundaries to the worldsheet surfaces. With this addition, we would like to perform the same topological twisting. To do so, we must first ensure that some supersymmetry condition persists when we put the boundary conditions on the fields. If we begin with a Calabi–Yau target manifold, we can ask to preserve the conditions that allow either the A-twisting or the B-twisting (but not both: the boundary condition will destroy some symmetry, much as pinning a rope will constrain its degrees of freedom). After the twist, the boundary topological theory will depend on symplectic or complex information, respectively.

For the A-model, the endpoints or boundaries must lie on a Lagrangian subspace. The Lagrangian condition constrains half the coordinates; for linear spaces it is like a restriction to the real part of a complex vector space. For the B-model the boundaries must lie on a complex space. Locally, a complex space looks like $\mathbb { C } ^ { n }$ and a complex subspace is described by complex analytic equations in the coordinates. A boundary condition that preserves supersymmetry and allows a chosen topological twisting is called a brane. (The terminology mimics the word “membrane,” but applies to any dimension.) In short, A-branes are Lagrangian; B-branes are complex.

To package all the information of the topological boundary theory, one appeals to the mathematical notion of a category [III.8]. A category is a way of talking about structure: it consists of objects, and for any pair of objects there is a space of morphisms from one object to the other. Often the objects are mathematical structures of some kind and the morphisms from one object to another are the functions that preserve the relevant structure. For example, if the objects are (i) sets [I.3 §2.1], (ii) topological spaces [III.90], (iii) groups [I.3 §2.1], (iv) vector spaces [I.3 §2.3], or (v) chain complexes, then the morphisms are, respectively, (i) maps [I.2 §2.2], (ii) continuous maps [III.90], (iii) homomorphisms [I.3 §4.1], (iv) linear maps [I.3 §4.2], or (v) chain maps. The morphism spaces between objects should be thought of as some kind of relational data. Morphisms themselves interact with one another, as they can be composed when the end object of one morphism is the start object of another. The composition is associative, so whether you compute abc as (ab)c or a(bc) does not matter. A useful image is a directed graph, which is a category with vertices as objects and paths between two vertices as morphisms. Composition is defined in this category by concatenating paths.

In the case of a two-dimensional field theory with boundary conditions, we construct a category whose objects are branes $( \mathrm { i . e . , }$ boundary conditions). The morphisms between two branes α and $\beta$ are the ground states $\mathcal { H } _ { \alpha \beta }$ of the boundary field theory defined on the infinite strip $[ 0 , 1 ] \times \mathbb { R } ,$ where we put the boundary condition α on the left boundary {0}×R and the condition β on the right boundary $\{ 1 \} \times \mathbb { R } .$ Morphisms are composed by gluing boundaries together, and associativity is guaranteed by topological invariance.14

Mirror symmetry with boundary conditions then becomes the following statement: two manifolds M and W are mirror pairs if the brane category of the Atwisting of M is equivalent to the brane category of the B-twisting of W (and vice versa). The mathematical translation of this statement is called the homological mirror symmetry conjecture, due to Kontsevich. On the A-model side, the brane category is the so-called Fukaya category, and is governed by complex analytic maps from surfaces with boundaries, where the boundaries must be mapped to Lagrangian branes. On the B-model side, the branes form a category determined by complex subspaces, together with complex analytic vector bundles [IV.6 §5] on them. A complex vector bundle associates a complex vector space to every point. For example, the complex circle $\{ x ^ { 2 } + y ^ { 2 } = 1 \}$ in $\mathbb { C } ^ { 2 }$ has a complex tangent space at every point. “Complex analytic” means that this subspace of $\mathbb { C } ^ { 2 }$ changes in a complex analytic way. For the complex circle, the space of tangent vectors at $( x , y )$ consists of all multiples of the vector $( - y , x )$ , an assignment which is clearly complex analytic. Physically, the bundles arise from allowing charges on the endpoints of strings.

Kontsevich’s conjecture asserts that these two categories of branes are equivalent. That statement is natural from the physics point of view, but by identifying the precise categories that correspond to the physical picture, this conjecture is a major contribution to the translation of mirror symmetry from physics into rigorous mathematics. The equivalence of categories means that not only is there a corresponding Lagrangian A-brane of M for every complex B-brane of W , but that the relationships, or morphisms, between branes are also in correspondence.

# 8.1 Example: Torus

Kontsevich’s conjecture can be proven and easily illustrated in the example of a 2-torus. Think of the nowfamiliar symplectic two-torus as being the two-dimensional plane, with integer lattice translations identified. We take the torus to have area element A dx dy, so that the symplectic parameter is the imaginary number $\rho = \mathrm { i } A$ , as in section 4.1.2. Now consider straight lines on the plane. These will correspond to closed circles on the torus as long as they have rational slope: $m = d / r$ , with d and r relatively prime integers. They are Lagrangian branes of the A-model boundary theory. The minimal-energy open strings connecting one line of slope m $= d / r$ to another of slope $m ^ { \prime } = d ^ { \prime } / r ^ { \prime }$ are those that have zero length. They are therefore the points of intersection. It is an easy exercise to show that there are $| d r ^ { \prime } - r d ^ { \prime } |$ such points.

On the mirror side, we again have a torus, but with a complex parameter τ, and for the two tori to be mirror pairs, we should set $\tau \ : = \rho$ . The objects of the Bmodel brane category are complex vector bundles. It is a theorem that the basic bundles are classified by their rank r and degree d, two integers.15 It is customary to organize these two numbers into what is known as a “slope,” $m \ = \ d / r$ (the nomenclature preceded this application), and basic bundles must have d and r relatively prime.

We can now easily guess that under the mirror correspondence we have

$\mathbf { s l o p e \longleftrightarrow { s l o p e } . }$

This means that a Lagrangian brane of slope m on the torus with symplectic parameter ρ should correspond to a complex vector bundle with slope m in the mirror torus with complex parameter $\rho .$ Now suppose we have the B-model version of our example above, so we take two vector bundles of slope m and $m ^ { \prime }$ . In fact, the minimum-energy open strings between two complex analytic bundles of slope m and $m ^ { \prime }$ correspond to complex maps between the bundles, and the riemann– roch formula [V.31] counts this number as $\ : | d \boldsymbol { r } ^ { \prime } - \boldsymbol { \mathbf { \ell } } | \ :$ r d	 . This is the same result as for our A-model calculation above! Therefore, corresponding objects relate in a corresponding way. Beyond the morphism spaces, one checks finally that the compositions of corresponding morphisms correspond, just as for logarithms and slide rules. Doing so proves Kontsevich’s conjecture.

# 8.2 Definition and Conjecture

In fact, Kontsevich’s definition of mirror symmetry is really a conjecture stating that the boundary notion of mirror symmetry as an equivalence of categories is compatible with, and even implies, the traditional notion of mirror symmetry that relates Gromov–Witten theory and complex structures.

One way to show this is to try to reconstruct the Gromov–Witten invariants from the boundary theory. A heuristic, geometric approach to doing so involves looking at the diagonal boundary condition in two copies of a space. A disk mapping into two copies of a space is described by two maps of a disk into the space. Further, if the boundary condition is diagonal, this means that the maps have to agree on the boundary. What we have, then, is two disks inside a space which agree on the boundary. That is exactly what a sphere is: two disks (or cups) glued together! The disks are the two hemispheres, and they are glued along the equator. Now the minimal disks are instantons for the open string (with boundary), and by gluing them together along a common boundary, we have constructed a minimal sphere, or closed-string instanton. Thus the open string on this double theory should recover the closed string on the original theory.

A more algebraic approach sees the closed-string deformations as deformations of the category of branes. That is, a change in bulk (nonboundary) theory induces a change in boundary theory. But once equipped with a category, one can classify its deformations intrinsically. That is, if one views a category as a fancy algebra,16 then, as the deformations of an algebra are easily classified through a notion called Hochschild cohomology, the deformations of a category can be treated similarly. One arrives at the maxim that the closed string is the Hochschild cohomology of the open string. By computing the Hochschild cohomology of a brane category, one can, in principle, check this maxim, establish Kontsevich’s conjecture, and then prove the connection to traditional mirror symmetry and Gromov–Witten theory.

# 9 Unifying Themes

How does one find mirror pairs (M, W)? What is the construction? Although mirror symmetry has spawned many results and proofs, these basic questions continue to vex.

On the one hand, Hori and Vafa have given a physics proof of mirror symmetry, which constructs mirror pairs but not through an evident mathematical channel. Of course, one can attempt to mathematicize the physical argument, but that does not seem to lead to insights into the construction—perhaps because path integrals and other methods of quantum field theory such as renormalization are not very well understood mathematically.

Batyrev has devised a procedure for constructing mirror pairs within the context of toric geometry. This method is a generalization, to a wide class of examples, of the original construction of Greene and Plesser. The recipe has been extremely successful in producing examples of every stripe. However, the underlying meaning behind the construction is unclear.

As for a geometric construction of mirror pairs, there is a physical argument that makes contact with mathematics, but it has not yet been made rigorous. The argument uses T-duality. Start with the B-model on M and consider a point P of M as a zero-dimensional complex subspace. Then the choice of point P on M is parametrized by M itself. By mirror symmetry, there should be a corresponding Lagrangian brane T on the mirror manifold W . Furthermore, the choices of T must equal the choices of P, i.e., the manifold M. Therefore, if we can find the brane T on W , we can parametrize the choices of T , and recover M. So we can find the mirror M of W from W itself.

This construction is geometric and has something to say about the structure of the Calabi–Yau spaces involved in mirror symmetry. Specifically, the choices of a Lagrangian brane always look like a family of tori. Therefore, M itself should look like a family of tori. Further, one can argue that by performing T-duality in families of tori (in a similar way to how one does it for a single torus), one arrives back at the mirror manifold, W . This is what we did for the torus, thought of as a circle family $( S _ { R _ { 2 } } ^ { 1 } )$ ) of circles $S _ { R _ { 1 } } ^ { 1 }$ . When we T-dualized each member of the family, we found the mirror torus. So mirror symmetry is T-duality, and Calabi–Yau spaces of mirror symmetry should look like families of tori. This approach also relates to the homological mirror symmetry construction. Though promising, it remains mathematically elusive.

Various points of view on mirror symmetry are helpful for different applications. To date, no unified understanding of the phenomenon has been achieved. To some extent, we are still “feeling the elephant.”

# 10 Applications to Physics and Mathematics

As a computational tool in string theory, mirror symmetry is unparalleled in its power. When combined with other physical equivalences, its power is multiplied. For example, there are certain equivalences in physics that relate one type of string theory to another.

Without going into the details of string theory, we can get a flavor of its complexity by returning to mirror symmetry. Recall that the B-model was able to compute the difficult instantons on the A-model, yielding a great simplification of the two-dimensional quantum field theory on the worldsheet. But this whole quantum field theory was just an auxiliary tool for computing some Feynman diagram for the perturbation theory of the full string theory! Unfortunately, a satisfactory description of the full string theory path integral is, at the time of writing, way out of reach. String theory instanton effects are mostly unknown to us, unless a string equivalence or other argument can relate them to a perturbative effect in a different string theory. The perturbative string calculation in that other theory may then be performed by exploiting mirror symmetry. Tracing through chains of equivalences in such a manner, many different phenomena in string theory can ultimately be calculated via mirror symmetry.

In principle, one should be able to calculate all nonperturbative and perturbative aspects of a single theory by outsourcing the calculations to equivalent theories and exploiting mirror symmetry. The barriers to doing this at the time of writing are largely technological, not conceptual.

Beyond physics, the rich texture of mirror symmetry means that there is interesting mathematics to be discovered in the proper formulation of the problem. For example, defining the precise categories of branes in full generality remains a challenge.

Yet there are also direct applications to mathematical questions. We have already discussed how enumerative geometry has been revolutionized by mirror symmetry and the counting of instantons. Results in symplectic geometry have also been obtained. Occasionally, two objects may be proven to be equivalent as B-model branes. If the A-model mirrors can then be found, one has the result that the corresponding Lagrangian subspaces of the mirror symplectic space are also equivalent. Of course, to make such an argument, one must first prove Kontsevich’s version of mirror symmetry for the mirror pair considered. As a final recent example, Kapustin and Witten have found a relation of mirror symmetry to the geometric Langlands program in representation theory. This program, loosely stated, is a correspondence between objects associated with twodimensional surfaces and Lie groups. From a surface Σ and a gauge group $G ,$ one constructs the space $\mathcal { M } _ { \mathrm { H } }$ of solutions to Hitchin’s equations. Central to that program are complex analytic objects on $\mathcal { M } _ { \mathrm { H } }$ that behave nicely under the action of an algebra of operations. The Langlands correspondence relates two sets of such objects: one easy to calculate and the other more difficult. In fact $\mathcal { M } _ { \mathrm { H } }$ is itself a family of tori, and the easy objects correspond to points. Mirror symmetry states that the points should turn into the tori under T-duality, so the hard objects should correspond to the tori themselves! It is an appealing proposition, and making it precise mathematics will be difficult—but the gauntlet has been thrown down.

The discovery that mirror symmetry relates to the geometric Langlands program has elicited great excitement among researchers and reveals yet another facet of this fascinating phenomenon.

# Further Reading

The article “Physmatics” (which can be found online at www.claymath.org/library/senior\_scholars/zaslow\_ physmatics.pdf) is a general discussion of the relationship between mathematics and physics, and may serve as a complement to this article. Readers with a university-level mathematics background who want to learn about mirror symmetry in more detail could try consulting the book Mirror Symmetry (Clay Mathematics Monographs, volume 1, edited by K. Hori and others (American Mathematical Society, Providence, RI, 2003)).

# IV.17 Vertex Operator Algebras Terry Gannon

# 1 Introduction

Algebra is the mathematics that places more emphasis on abstract structure than on intrinsic meaning. The conceptual simplifications that can result when context is stripped away from structure give algebra a special power and clarity compared with other areas: compare, for example, the difficulty of visualizing fourdimensional space with the triviality of manipulating quadruples $( x _ { 1 } , x _ { 2 } , x _ { 3 } , x _ { 4 } )$ of real numbers. However, this abstractness can also blind us. For instance, basic identities like ab  ba and $\begin{array} { l l l } { { a ( b c ) } } & { { = } } & { { ( a b ) c } } \end{array}$ that are obeyed by numbers can be modified in countless directions, and each modification defines a new algebraic structure, but it is hard to guess from a purely abstract perspective which of these modifications will give rise to a rich, accessible, and interesting theory. For guidance, algebra has traditionally turned to geometry. For example, over a century ago lie [VI.53] suggested that the identities ab = −ba and a(bc) = $( a b ) c + b ( a c )$ were worth studying for geometrical reasons: the resulting structures are now called lie algebras [III.48 §2]. More recently, as we shall see, physics has joined geometry in this guiding role and has had spectacular success.

The renowned physicist and mathematician Edward Witten believes that a major theme of twenty-firstcentury mathematics will be its reconciliation with the branch of physics known as quantum field theory. Conformal field theory (the quantum field theory that underlies string theory) is an especially symmetric and well-behaved class of quantum field theories. When this notion is translated into algebra, the result is a structure known as a vertex operator algebra (VOA). This article sketches where VOAs come from, what they are, and what they are good for.

To aim to explain a VOA in a few pages is almost as absurd as to aim to explain quantum field theory in a few pages, but, undaunted, I shall try to do both. Obviously it will be necessary to gloss over many important technicalities and to commit major simplifications; without question this exposition will raise the ire of experts and the eyebrows of knowledgeable amateurs, but I hope that it will at least convey the essence of this important and beautiful area. Vertex operator algebras are the algebra of string theory: they should be thought of as the same sort of gift to the twenty-first century that Lie algebras were to the twentieth.

# 2 Where VOAs Come From

The two most revolutionary developments in physics in the early twentieth century are usually held to be relativity and quantum mechanics. They are revolutionary not just because they have consequences that are extremely counterintuitive, but also because they provide very general frameworks that can potentially affect all physical theories: one can take a theory from classical physics, such as the theory of the harmonic oscillator or the theory of electrostatic force, for example, and one can try to make it “relativistic,” so that it becomes compatible with relativity, or to “quantize” it, so that it becomes compatible with quantum mechanics.

Unfortunately, nobody knows how to make relativity fully compatible with quantum mechanics. To put this another way, the ultimate concern of relativity is gravitation, and a direct application to gravity of the usual quantizing techniques fails. This ought to mean that a fundamentally new physics arises at small distance scales that we are ignoring. Indeed, naive calculations suggest that the space-time “continuum” at distance scales of around $1 0 ^ { - 3 5 }$ m should deteriorate into some sort of “quantum foam,” whatever that might mean. $( 1 0 ^ { - 3 5 }$ m is extremely small: for instance, the order of magnitude of the size of an atom is $1 0 ^ { - 1 0 }$ m.)

Perhaps the most popular and controversial approach to quantum gravity is string theory. The electron is a particle, i.e., in principle it can be localized to a point. In string theory, the fundamental object is a string, a finite curve of length approximately $1 0 ^ { - 3 5 }$ m. In place of the dozens of kinds of fundamental particles in the generally accepted quantum field theory, there is only one string, whose precise physical properties (mass, charge, etc.) depend on its current “vibrational mode.”

As the string moves, it traces out a surface called a worldsheet. For reasons that we will sketch below, much of string theory reduces to studying conformal field theory, which is the induced quantum theory on these surfaces. Probably no other structures have affected so many areas of “pure” mathematics in so short a time as string theory and, what is essentially the same thing, conformal field theory. Indeed, five of the twelve Fields Medals awarded in the 1990s (namely, those to Drinfel’d, Jones, Witten, Borcherds, and Kontsevich) were for such work. We shall focus in this article on their algebraic impact; see mirror symmetry [IV.16] for some geometrical implications.

# 2.1 Physics 101

A quick overview of physics will be useful for the discussion. Further details can be found in mirror symmetry [IV.16 §2].

# 2.1.1 States, Observables, and Symmetries

A physical theory is a set of laws that govern the behavior of some kind of physical system. A state of that system is a complete mathematical description of the system at a particular time: for instance, if the system consists of a single particle, then we could take its state to be its position x and momentum $p \ = \ m ( { \bf d } / { \bf d } t ) x$ (where m is its mass). An observable is a physically measurable quantity such as position, momentum, or energy. It is through observables that a theory is compared with experiment. Of course, for this to be true we also need to know what an observable is from a theoretical point of view.

In classical physics, an observable is just a numerical function of the state. For example, our single particle has energy E, which depends on the position and momentum via a formula of the form E $( 1 / 2 m ) p ^ { 2 } + V ( x )$ . (This gives us the kinetic energy plus the potential energy.) Classical states at different times are related by the equations of motion, which are usually expressed as differential equations. However, string theory and conformal field theory (CFT) are quantum theories, which are significantly different from classical theories: one can think of them as “applied linear algebra.” Whereas a classical state was given by a collection of a few numbers (two, in the case of the particle above), a quantum state is an element of a hilbert space [III.37], which for the purposes of discussion we can think of as a column vector with infinitely many complex entries. As for a quantum observable, it is a hermitian operator [III.50 §3.2] on the Hilbert space, which we can think of as an  matrix Aˆ that acts on the states by matrix multiplication. As in classical physics, one of the most important observables is energy, which is given by the Hamiltonian operator Hˆ.

It is far from obvious how a linear operator that takes states to states has anything to do with the notion of a physical observation, and indeed the relationship between observables and observation is a major difference between classical and quantum theories. If Aˆ is an observable, then the spectral theorem [III.50 §3.4] tells us that the Hilbert space has an orthonormal basis [III.37] of eigenvectors [I.3 §4.3]. When we do the experiment that is modeled by the observable Aˆ, the answer we obtain will be one of the eigenvalues of Aˆ. However, this answer is usually not fully determined by the state v. Instead, it is given by a probability distribution: the probability of obtaining a particular eigenvalue is proportional to the square of the norm of the projection of v into the corresponding eigenspace. Thus, the only circumstances under which the answer is determined in advance are if the state v is an eigenvector of Aˆ.

There are two independent ways in which a quantum state can evolve in time: a deterministic evolution between measurements, governed by the famous schrödinger equation [III.83], and a probabilistic and discontinuous one that occurs at the instant when a measurement is made. For our purposes, only the deterministic evolution will be relevant.

The symmetries of CFT are extremely rich, as we shall see. Symmetries in physical theories are highly desirable because of two consequences that they have. First, they lead by noether’s theorem [IV.12 §4.1] to conserved quantities, i.e., quantities independent of time. For example, the equations of motion of our particles are usually invariant under translation: for instance, the gravitational force between two particles depends only on the difference between their positions. The corresponding conservation law in this case is the conservation of momentum. A second consequence of symmetries in quantum theories is that infinitesimal generators of the symmetries act on the state space H (the Hilbert space to which the states belong), forming a representation of the Lie algebra. Both consequences are important to CFT.

# 2.1.2 The Lagrangian Formulation and Feynman Diagrams

We will need two of the languages in which physics is written. One is the Lagrangian formalism, which is responsible for the relationship between string theory and CFT, as well as for the appearance of modular functions in string theory. The other is the Hamiltonian or Poisson bracket formalism, which is where algebra arises. Vertex operator algebras try to explain the “miracle” that these two formalisms cohere.

The Lagrangian formalism can be expressed classically through Hamilton’s action principle. When there are no forces present, particles travel in straight lines, which are the curves of shortest length. Hamilton’s principle explains how this idea generalizes to arbitrary forces: instead of minimizing length, the particle minimizes a related quantity S called the action.

The quantum version of Hamilton’s principle is due to Feynman. He expresses the probability of measuring the system in some final (eigen)state |out, given that it was originally in some initial state |in, using a “path integral” of eiS/ over all possible histories that connect in and out . The details are not important for us (and in any case are mathematically dubious in general). The intuition behind the path integral formulation is that the particle simultaneously follows every one of those histories, and each of them contributes to the probability.  is called Planck’s constant; in the “classical limit” as   0, the contribution from the path that satisfies Hamilton’s principle dominates everything else.

The main use of Feynman’s path integral is in perturbation theory. Finding exact solutions in physics is typically impossible and rarely useful. In practice, it suffices to find the first few terms in some Taylor expansion of the solution. This so-called “perturbative” approach to quantum theories is particularly transparent in Feynman’s formalism, where each term of the expansion can be represented pictorially as a graph. See figure 1(a) for typical examples. The graphs contributing to the nth-order term in this Taylor expansion will involve n vertices. Feynman’s rules describe how to convert these graphs into integral expressions for computing the individual terms in the Taylor expansion.

In this article we are interested in perturbative string theory. The string Feynman diagrams (see figure 1(b) for three equivalent ones) are surfaces called worldsheets; the need for quantum foam is avoided because these surfaces are much less singular than the particle graphs (which have singularities at each vertex), and this is also largely why the mathematics of strings is so nice. To cut a long story short, each term in the perturbative expression for probabilities in string theory can be calculated from a quantity called a “correlation function” in a CFT that lives on the corresponding

![](images/eed7f334947b0752686c514e1581c1faf9fa2c3f2d49e3415507c11d41fa709c.jpg)

<details>
<summary>natural_image</summary>

Diagram showing two abstract geometric structures labeled (a) and (b), with no text or symbols present.
</details>

Figure 1 Some Feynman diagrams of (a) particles and (b) strings.

worldsheet. Feynman’s path integral here amounts to the integral of a quantity that CFT can compute, over some moduli space [IV.8] of surfaces.

The vertices in a Feynman diagram represent places where one particle absorbs or emits another. The corresponding rules of string theory tell us that we should dissect the worldsheet into “tubular Y-shapes,” or spheres with three legs, as in figure 2. Since these spheres with legs play the role of vertices in the Feynman diagram, the factor they contribute to the integrand of the path integral is called a vertex operator, and now it describes the absorption or emission of one string by another. A vertex operator algebra is the “algebra” of these vertex operators.

# 2.1.3 The Hamiltonian Formulation and Algebra

The Poisson bracket A, B P of two classical observables A and B is defined to be

$$
\frac {\partial A}{\partial x} \frac {\partial B}{\partial p} - \frac {\partial B}{\partial x} \frac {\partial A}{\partial p}.
$$

Note that $\{ A , B \} _ { \mathrm { P } } = - \{ B , A \} _ { \mathrm { P } } \{$ in other words, the Poisson bracket is anti-commutative. It also satisfies the Jacobi identity

$$
\{A, \{B, C \} _ {\mathrm{P}} \} _ {\mathrm{P}} + \{B, \{C, A \} _ {\mathrm{P}} \} _ {\mathrm{P}} + \{C, \{A, B \} _ {\mathrm{P}} \} _ {\mathrm{P}} = 0,
$$

and therefore defines a Lie algebra. The Hamiltonian formulation of classical physics expresses the evolution of an observable A by means of the differential equation $\dot { A } \ = \ \{ A , H \} _ { \mathrm { { P } } }$ , where H is the hamiltonian [III.35]: that is, the energy observable. The quantum version of this picture is due to Heisenberg and Dirac: the observables are now linear operators rather than smooth functions, and the Poisson bracket is replaced by the commutator $[ \hat { A } , \hat { B } ] = \hat { A } \circ \hat { B } - \hat { B } \circ \hat { A }$ of operators. This again has the anti-commuting property $[ \hat { A } , \hat { B } ] =$ $[ \hat { B } , \hat { A } ]$ and again satisfies the Jacobi identity, so the process of “quantization” gives rise to a homomorphism of Lie algebras. The derivative with respect to time of a quantum observable Aˆ is then the natural analogue of the classical case: it is proportional to $[ \hat { A } , \hat { H } ]$ , where Hˆ is the Hamiltonian operator. Thus the Hamiltonian has a dual role: as the energy observable and as the controller of time evolution. All of physics is stored in the action of the observables on state space ${ \mathcal { H } } ,$ as well as the commutators of these observables with Hˆ.

![](images/59c0018bd1f92a6ebd6e8d542e8f8c88735dfa1e4eaaeef37e33bceedd6ec2e0.jpg)

<details>
<summary>natural_image</summary>

Two abstract line drawings of symmetrical, organic shapes with internal cutouts (no text or symbols)
</details>

Figure 2 Dissecting a surface.

Let us illustrate this picture with the quantum spring, also known as the harmonic oscillator. The position and momentum observables $\hat { x } , \hat { p }$ are operators acting on the infinite-dimensional space  of possible springstates. It is more convenient to work with certain combinations of them called aˆ and $\hat { a } ^ { \dagger }$ (the dagger denotes the “Hermitian adjoint,” or complex-conjugate transpose), which obey the commutator relation $[ \hat { a } , \hat { a } ^ { \dagger } ] = I ,$ , where I is the identity operator. It turns out that all other observables can be built from aˆ and $\hat { a } ^ { \dagger }$ . For example, the Hamiltonian Hˆ is $l ( \hat { a } ^ { \dag } \hat { a } + \frac { 1 } { 2 } )$ for some positive constant l. The vacuum, which is denoted |0, is the state of minimum energy. In other words, the state 0 is an eigenvector of $\hat { H }$ with smallest possible eigenvalue: $\hat { H } | 0 \rangle = E _ { 0 } | 0 \rangle$ for some $E _ { 0 } \in$ R and all other eigenvalue E of Hˆ are greater than $E _ { 0 } .$ . It follows from this that ${ \hat { a } } | 0 \rangle = 0 .$ . To see why, consider the effect of $\hat { H }$ on aˆ|0:

$$
\begin{array}{l} \hat {H} \hat {a} | 0 \rangle = l (\hat {a} ^ {\dagger} \hat {a} + \frac {1}{2}) \hat {a} | 0 \rangle = l (\hat {a} \hat {a} ^ {\dagger} - \frac {1}{2}) \hat {a} | 0 \rangle \\ = \hat {a} l (\hat {a} ^ {\dagger} \hat {a} - \frac {1}{2}) | 0 \rangle = \hat {a} (\hat {H} - l) | 0 \rangle = (E _ {0} - l) \hat {a} | 0 \rangle . \\ \end{array}
$$

Here, we have used the fact that $\hat { a } ^ { \dagger } \hat { a } = \hat { a } \hat { a } ^ { \dagger } - I .$ (The observables aˆ and $\hat { a } ^ { \dagger }$ are called creation and annihilation operators because, as we shall see later, they can be interpreted as adding or removing a particle from a certain n-particle state. Showing this uses the fact that they produce I when you interchange their order.) This calculation shows that if $\hat { a } | 0 \rangle$ is not zero, then it is an eigenvector of Hˆ with an eigenvalue smaller than $E _ { 0 } ,$ , which is a contradiction.

Since ${ \hat { a } } | 0 \rangle = 0 ,$ it follows that $\hat { H } | 0 \rangle = \textstyle \frac { 1 } { 2 } l | 0 \rangle$ , so $E _ { 0 } =$ ${ \frac { 1 } { 2 } } l .$ We now define, for each positive integer n, a state |n to be $( \hat { a } ^ { \dag } ) ^ { n } | 0 \rangle \in \mathcal { H }$ . Similar calculations to the one just given show that |n has energy $E _ { n } = ( 2 n + 1 ) E _ { 0 }$ . For example,

$$
\begin{array}{l} \hat {H} | 1 \rangle = l (\hat {a} ^ {\dagger} \hat {a} + \frac {1}{2}) \hat {a} ^ {\dagger} | 0 \rangle = l (\hat {a} ^ {\dagger} (\hat {a} ^ {\dagger} \hat {a} + I) + \frac {1}{2} \hat {a} ^ {\dagger}) | 0 \rangle \\ = \frac {3}{2} l \hat {a} ^ {\dagger} | 0 \rangle = E _ {1} | 1 \rangle . \\ \end{array}
$$

(Note that we used the fact that $a | 0 \rangle = 0$ in the penultimate equality above.) We think of the vacuum as the ground state, and n as being the state with n quantum particles. These states n span all of the state space . To see how some observable acts on some state, one writes the observable in terms of the basic observables $\hat { a } , \hat { a } ^ { \dag }$ and the state in terms of the basic states n . In this algebraic way we can recover all of the physics.

This idea of building up the whole space from the vacuum and the operators is a fruitful one in mathematics as well: something similar happens for the most important modules of most of the important Lie algebras.

# 2.1.4 Fields

A classical field is a function of space and time. Its values can be numbers or vectors, which represent quantities such as air temperature or the current in a river. The values taken by a quantum field are operators; furthermore, a quantum field is not a function of space and time, but a more general object called a distribution [III.18]. The prototypical example of a distribution is the Dirac delta function $\delta ( x - a )$ . Despite its name, this is not a function: rather, it is defined by the property that

$$
\int f (x) \delta (x - a) \mathrm{d} x = f (a) \tag {1}
$$

for any sufficiently well-behaved function $f ( x )$ . Even though $\delta ( x - a )$ is not a function, one can informally interpret it as the derivative of a step function, and one can visualize it as equaling 0 everywhere except at $x \ = \ a$ , where it is infinite, in such a way that the infinitely tall and infinitely thin rectangle under the graph has area 1. However, it really only makes sense inside an integral, as in (1). Similar remarks apply to distributions in general, so a quantum field can really only be evaluated inside an integral of space and time, applied to some “test function” like f above. The value of such an integral will be an operator on the state space .

Dirac deltas appear in classical mechanics when one takes Poisson brackets of classical fields. Similarly, commutators of quantum fields involve delta functions too. For example, in the simplest cases the quantum fields ϕ satisfy

$$
\left. \begin{array}{c} [ \varphi (x, t), \varphi (x ^ {\prime}, t) ] = 0, \\ \left[ \varphi (x, t), \frac {\partial}{\partial t} \varphi (x ^ {\prime}, t) \right] = \mathrm{i} \hbar \delta (x - x ^ {\prime}). \end{array} \right\} \tag {2}
$$

This is a mathematical way of expressing, in the context of quantum field theory, the cherished physical principle called locality: 1 the only way we can directly affect something is by nudging it. In order to influence something not touching us, we must propagate a disturbance from us to it, such as a ripple in water. The main purpose of both classical and quantum fields is that they provide a natural vehicle for realizing locality. Locality is also at the heart of vertex operator algebras.

An important aspect of modern physics is that many of the central concepts of classical physics become less central, and are instead derived quantities. For example, the basic object of general relativity [IV.13] is a Lorentzian manifold, and familiar physical quantities such as mass and gravitational force are, from the point of view of this manifold, just names (that are not wholly precise) given to certain of its geometrical features.

Particles are obviously essential to classical physics, but we have not mentioned them in our brief sketch of quantum field theory. They arise through the so-called modes of quantum fields ϕ, which play the role of the operators aˆ, aˆ† that we met in section 2.1.3. A mode is the operator that results from hitting the quantum field with an appropriate test function and integrating—just as one does when working out a Fourier coefficient, in which case the test functions are trigonometric functions [III.92]. In fact, when viewed appropriately, modes actually are Fourier coefficients of a certain kind. The commutators of these modes can be obtained from the commutators of the fields. Now, recall that the vertex operators of string theory are related to the emission and absorption of strings. As we shall see shortly, these vertex operators are the quantum fields in a quantum field theory of point particles (namely, the associated conformal field theory); the modes of these vertex operators generate the “particles” (or in more conventional language, the states) in that conformal field theory. Equivalently, they generate the various vibrational states of a single string in that string theory.

# 2.2 Conformal Field Theory

A conformal field theory (CFT) is a quantum field theory with a two-dimensional space-time whose symmetries include all conformal transformations. We shall explain what this means in the next paragraphs, but for now it is enough to know that a CFT is a particularly symmetrical kind of quantum field theory. A CFT lives on the worldsheet Σ traced by a set of strings as they evolve, sometimes colliding and separating, through time. In this subsection we shall informally sketch their basic theory; in section 3.1 we shall be more precise.

CFT, like any quantum field theory in two dimensions, has two almost independent halves. This is easiest to see in the context of string theory: the ripples on the string are responsible for the physical properties (charge, mass, etc.) of the corresponding state, but they can move (at the speed of light) either clockwise or counterclockwise around the string. When they do so, they just pass through each other without interacting. These two alternatives, clockwise and counterclockwise, yield the two chiral halves of CFT. To study a CFT, one first analyzes its chiral halves and then splices them together to form the “bichiral” physical quantities. Almost all attention in CFT by mathematicians has focused on the chiral (as opposed to physical) data, and indeed that is where vertex operator algebras live. For ease of presentation, we will usually suppress one of the chiral halves.

A conformal transformation is a transformation that preserves angles. The simplest reason one can give for why two dimensions are so special for CFT is that there are far more conformal transformations in two dimensions than there are in higher dimensions. When n > 2 the only examples are the obvious ones: combinations of translations, rotations, and enlargements. This means that the space of all local conformal transformations in $\mathbb { R } ^ { n }$ is $\binom { \binom { 1 } { n + 2 } } { 2 }$ dimensional. However, when n 2 the space of local conformal transformations is far richer: it is infinite dimensional. Indeed, if you identify $\mathbb { R } ^ { 2 }$ with the complex plane C, then any holomorphic function [I.3 §5.6] f (z) that does not have zero derivative at a point $z _ { 0 }$ is conformal near $z _ { 0 }$ . Since a CFT is invariant under conformal transformations and there are many conformal transformations, a CFT is especially symmetrical: this is what makes CFTs so interesting mathematically.

Lie algebras arise naturally whenever one has local symmetries, and indeed one can form an infinitedimensional Lie algebra out of the infinitesimal conformal transformations. This algebra has a basis $l _ { n } , n \in \mathbb { Z } ,$ that obeys the Lie-bracket relations

$$
[ l _ {m}, l _ {n} ] = (m - n) l _ {m + n}. \tag {3}
$$

The algebraic interpretation of the conformal symmetry of CFT turns out to be that these basis elements $l _ { n }$ act naturally on all the quantities in the theory, as we shall explain below.

The basic example that underlies all the others is when space-time Σ is a semi-infinite cylinder corresponding to an incoming string. It is parametrized by time $t < 0$ and the angle $0 \leqslant \theta < 2 \pi$ around the string. We can conformally map the cylinder to the punctured disk in C by $z ~ = ~ \mathrm { e } ^ { t - \mathrm { i } \theta }$ , so $t ~ = ~ - \infty$ corresponds to $z = 0 .$ . This allows us to say what we mean by conformal symmetries of the cylinder.

The quantum fields $\varphi ( z )$ of CFT are the vertex operators of string theory. As always, these quantum fields ϕ are “operator-valued distributions” on space-time $\Sigma ,$ acting on the space  of states. Now it is possible for a field ϕ to be “holomorphic,” in the following sense. First, you calculate its modes $\varphi _ { n }$ , one for each $n \in \mathbb { Z } ,$ , which are linear maps from the state space to itself, given by the formula

$$
\varphi_ {n} = \int \varphi (z) z ^ {n - 1} \mathrm{d} z,
$$

where the integral is around a small circle about the origin. Then you take these modes as the coefficients of a formal power series $\scriptstyle \sum _ { n \in \mathbb { Z } } \varphi _ { n } z ^ { n }$ . We call ϕ holomorphic if this formal power series can be identified with ϕ, in a sense that we shall discuss more in section 3.1. A typical field ϕ(z) is not holomorphic: rather, it is a combination of holomorphic and anti-holomorphic fields, which make up the two chiral halves of CFT. We will focus on the space of holomorphic fields ϕ(z), which we call  . This turns out to form a vertex operator algebra (as do the anti-holomorphic fields).

For example, the most important vertex operator comes directly from the conformal symmetry: the stress-energy tensor $T ( z ) \in \mathcal { V }$ is the “conserved current” that Noether’s theorem associates with the conformal symmetry. Labeling its modes (Noether’s “conserved charges” here) by $L _ { n } = \int T ( z ) z ^ { - n - 3 }$ dz, so that $\begin{array} { r } { T ( z ) = \sum _ { n } L _ { n } z ^ { - n - 2 } , } \end{array}$ , we find that they almost realize the conformal algebra: instead of (3), however, they obey the slightly more complicated relations

$$
[ L _ {m}, L _ {n} ] = (m - n) L _ {m + n} + \delta_ {n, - m} \frac {m (m ^ {2} - 1)}{1 2} c I, \tag {4}
$$

where I is the identity. In other words, the operators $L _ { n }$ and I form an extension of the conformal algebra by I. The resulting infinite-dimensional Lie algebra is called the Virasoro algebra . The number c appearing in Vir(4) is called the central charge of the CFT and is a rough measure of its size.

The operators $L _ { n }$ do not precisely represent the conformal algebra (3). Instead, they form a so-called projective representation. Projective representations of symmetries, such as (4), are common in quantum theories. The fact that they are not true representations is not a problem, since one can turn them into true representations by extending the algebra. In our case, the state space  carries inside it a true representation of the Virasoro algebra , which is useful as it means Vircan be used to organize .

Any quantum field theory has what is called a state– field correspondence: with each field ϕ one associates its incoming state, which is the limit as the time t tends to of ϕ 0 (as always, 0 is the vacuum state in and ϕ acts on states). CFT is unusual in that the state– field correspondence is a bijection. This means we can identify  and  and use states to label all fields.

We want to make V into some sort of algebra, but the obvious direct approach of taking products $\phi _ { 1 } ( z ) \phi _ { 2 } ( z )$ fails, since distributions, unlike true functions, cannot in general be multiplied. For example, the Dirac delta $\delta ( x - a )$ cannot be squared without causing problems in (1). However, even if the product $\phi _ { 1 } ( z ) \phi _ { 2 } ( z )$ does not make sense, one can make sense of $\varphi _ { 1 } ( z _ { 1 } ) \varphi _ { 2 } ( z _ { 2 } )$ as an operator-valued distribution on $\Sigma ^ { 2 }$ . It is then possible to recover most of the physics of CFT by studying the singular terms as $z _ { 2 } ~  ~ z _ { 1 }$ . By the operator product expansion, we mean expanding products $\varphi _ { 1 } ( z _ { 1 } ) \varphi _ { 2 } ( z _ { 2 } )$ as sums of the form $\Sigma _ { h } ( z _ { 1 } - z _ { 2 } ) ^ { h } O _ { h } ( z _ { 1 } )$ . The set is closed under this product in the sense that each coefficient $O _ { h } ( z )$ lies in  . A typical example is

$$
\begin{array}{l} T (z _ {1}) T (z _ {2}) = \frac {1}{2} c (z _ {1} - z _ {2}) ^ {- 4} I + 2 (z _ {1} - z _ {2}) ^ {- 2} T (z _ {1}) \\ + \left(z _ {1} - z _ {2}\right) \frac {\mathrm{d}}{\mathrm{d} z} T \left(z _ {1}\right) + \dots . \\ \end{array}
$$

Physicists call a chiral algebra; for us it is the prototypical example of a vertex operator algebra. It is not an algebra in the conventional sense though, since, given vertex operators $\varphi _ { 1 } ( z )$ and $\varphi _ { 2 } ( z )$ , we have not just a single product $\varphi _ { 1 } ( z ) * \varphi _ { 2 } ( z )$ in V but infinitely many products $\varphi _ { 1 } ( z ) * _ { h } \varphi _ { 2 } ( z ) = { O } _ { h } ( z )$ , all belonging to $\mathcal { V } .$ .

The Hamiltonian plays a crucial role in any quantum field theory; here it turns out to be proportional to the mode $L _ { \mathrm { 0 } }$ discussed earlier. Being an observable, $L _ { \mathrm { 0 } }$ is diagonalizable on ${ \mathcal { H } } ,$ , which means that any state $\nu \in$ H can be written as a sum $\sum _ { h } \nu _ { h }$ , where $\nu _ { h } \in \mathcal { H }$ has energy h: that is, $L _ { 0 } \nu _ { h } = h \nu _ { h }$ .

There is a special class of CFT that is particularly well-behaved. Let ¯ denote the space of all anti-holomorphic fields in the CFT—it is the other chiral half. Recall that the full CFT consists of  and ¯ spliced together. We call the CFT rational if $\mathcal { V } \oplus \bar { \mathcal { V } }$ is so large that it has finite index, in an appropriate sense, in the full space of quantum fields in the CFT. The name “rational” arises because the central charge c and other parameters in a rational CFT have to be rational numbers.

The mathematics of rational CFT is especially rich. Let us briefly look at one example. (We will use several words that will be unfamiliar to most readers, but at least it will give some idea of which areas are touched by CFT.) As with everything else, the quantum probabilities arising in CFT are found by first computing chiral quantities and splicing them together. These chiral quantities are called conformal or chiral blocks, and are found using simple Feynman-like rules applied to dissections like figure 2. In rational CFT we get a finitedimensional space ${ \mathcal { F } } _ { g , n }$ of chiral blocks for any worldsheet $\Sigma , \mathbf { i . e . }$ , for any choice of genus $_ g$ and number n of punctures. These spaces carry projective representations of the mapping class group $T _ { g , n }$ (defined to be the fundamental group $\pi _ { 1 }$ of the moduli space ${ \mathcal { M } } _ { g , n } ) .$ . This $T _ { g , n }$ -representation is the source, for instance, of Jones’s relation of the braid group [III.4] (and hence knots [III.44]) to subfactors, Borcherds’s explanation of “Monstrous Moonshine,” the Drinfel’d–Kohno monodromy theorem, and the modularity of affine Kac– Moody characters. Some of this we will touch on in section 4.

The most important example here is the torus, where the chiral blocks are modular functions, a class of functions of fundamental mathematical importance. A modular function is a meromorphic function (that is, a function that is holomorphic except at a few “poles” where it can tend to infinity) f (τ) that is defined on the upper half-plane $\mathbb { H } = \{ \tau \in \mathbb { C } \mid \operatorname { I m } \tau > 0 \}$ and that is “symmetric” with respect to the group SL2(Z) of $2 \times 2$ matrices with integer entries and determinant 1, in the sense that for any such matrix $( \mathbf { \Pi } _ { c } ^ { a \ b } )$ ) the function $f ( \tau )$ is closely related (though not necessarily exactly equal) to the function $f ( ( a \tau + b ) / ( c \tau + d ) )$ ). We shall discuss this further in section 3.2.

The appearance of modularity can be understood by recalling from section 2.1.2 that Feynman’s path integral in string theory is an integral over moduli spaces. The moduli space $\mathcal { M } _ { 1 , 0 }$ for the torus can be written as the quotient of the half-plane H by the action of $\mathrm { S L } _ { 2 } ( \mathbb { Z } )$ . Therefore, if one lifts the integrand of Feynman’s integral from $\mathcal { M } _ { 1 , 0 }$ to H, one obtains a function $\mathcal { Z } ( \tau )$ that is invariant under $\mathrm { S L } _ { 2 } ( \mathbb { Z } )$ and hence modular. This integrand (τ) is a quadratic combination of the chiral blocks for the torus.

# 3 What VOAs Are

It is possible to give a fully axiomatic definition of vertex operator algebras. However, when one first encounters this definition (and not just the first time either) it can seem very complicated and arbitrary, and one is given no feel for the importance of VOAs. Our treatment below will be much more informal: this will clarify their importance even if it hides much of their complexity. Thanks to the previous section, it is possible to give a quick justification for VOAs: if you concede that CFT (or equivalently, perturbative string theory) is important, and if you have seen how closely related CFT is to VOAs, then you must concede that VOAs are important. However, this is not the whole story, as we shall see.

# 3.1 Their Definition

Let us begin by defining them in terms of other concepts that must themselves be defined: a vertex operator algebra is an algebra of vertex operators, or in other words the chiral algebra V of a conformal field theory.

The most important thing to understand in this definition is that a vertex operator is a quantum field, which, as we have seen, is an “operator-valued distribution of space-time.” So we can think of it informally as a matrix-valued function of space-time, where the matrix is $\infty \times \infty$ and its entries can be generalized functions like the Dirac delta (1). However, we shall give a much better description of these vertex operators shortly.

By “space-time” we mean the unit disk in C punctured at $z = 0$ . Recall from section 2.2 that string-theoretically this set corresponds to a semi-infinite cylinder parametrized by the angle $- \pi < \theta \leqslant \pi$ running around the string as well as the time $\vdash \infty \ < \ t \ < \ 0$ running along the axis: the map from this to the punctured disk was $( \theta , t ) \mapsto z = \mathrm { e } ^ { t - \mathrm { i } \theta }$ . We want to restrict our attention to quantum fields that depend holomorphically on z. However, it is not obvious what “holomorphic” means for distributions. We touched on this question in section 2.2: now we shall look at it in more detail.

To do this, we need a more concrete description of a vertex operator. The key idea is a very convenient algebraic interpretation of holomorphic distributions. Consider the sum

$$
d (z) = \sum_ {n = - \infty} ^ {\infty} z ^ {n}. \tag {5}
$$

Multiply it by $f ( z ) = 3 z ^ { - 2 } - 5 z ^ { 3 }$ , say. This gives us

$$
\begin{array}{l} f (z) d (z) = 3 \sum_ {n = - \infty} ^ {\infty} z ^ {n - 2} - 5 \sum_ {n = - \infty} ^ {\infty} z ^ {n + 3} \\ = 3 \sum_ {n = - \infty} ^ {\infty} z ^ {n} - 5 \sum_ {n = - \infty} ^ {\infty} z ^ {n} = - 2 d (z). \\ \end{array}
$$

A few more examples like this will convince you that $f ( z ) d ( z ) = f ( 1 ) d ( z )$ for any polynomial function $f$ of z and $z ^ { - 1 }$ . Therefore, d(z) behaves exactly like the Dirac delta $\delta ( z - 1 )$ , at least for polynomial test functions f . Note that d(z) cannot converge for any z: the positive powers have a convergent sum only for $| z | < 1$ , and the negative powers only for $| z | > 1$ . The “function” d(z) is an example of a formal power series: any series $\scriptstyle \sum _ { n = - \infty } ^ { \infty } a _ { n } z ^ { n }$ , where the coefficients $a _ { n }$ can be anything and we ignore all convergence issues.

By inspection, these formal power series are “holomorphic” throughout the punctured plane: after all, holomorphic just means that the complex derivative d/dz exists, and the derivative $\scriptstyle \sum _ { n } n a _ { n } z ^ { n - 1 }$ of a formal power series clearly remains a formal power series. (By contrast, nonholomorphic series would involve the complex conjugate z¯.)

So that is what a vertex operator looks like: a formal power series $\scriptstyle \sum _ { n = - \infty } ^ { \infty } a _ { n } z ^ { n }$ , where each coefficient $a _ { n }$ is now an operator (endomorphism) on the space of states, which is an infinite-dimensional vector space. Since the vertex operators are in one-to-one correspondence with the states (we called this the “state–field correspondence” above), we can label these vertex operators with states: the standard convention is to denote the vertex operator corresponding to state $\upsilon \in \mathcal { V }$ by

$$
Y (\nu , z) = \sum_ {n = - \infty} ^ {\infty} \nu_ {n} z ^ {- n - 1}. \tag {6}
$$

The symbol $" Y "$ should remind you of the sphere with three legs, which as we know is the vertex of string theory. These coefficients $\nu _ { n }$ are the modes: as in any quantum field theory, all observables and states in the theory are built up from them.

The most important state in the theory is the vacuum |0. It corresponds to the identity vertex operator: $Y ( \left| 0 \right. , z ) = I .$ . From the physical point of view, the vertex operator $Y ( \nu , z )$ is the field that created the state v at time $t = - \infty$ , i.e., Y (v, 0) 0 exists and equals v. (Recall that in our model z  0 corresponds to $t = - \infty . )$ Among other things, this means that $\nu _ { - 1 } ( | 0 \rangle ) ~ = ~ \nu _ { \mathrm { { } } }$ , so indeed the modes applied to 0 generate $\gamma ,$ as is required in any quantum field theory.

The most important observable in the theory is the Hamiltonian, or energy operator, which we denote by $L _ { 0 }$ . It is diagonalizable (so V can be written as a sum of L0-eigenspaces) and all of its eigenvalues must be integers. For example, the vacuum |0 has 0 energy: $L _ { 0 } | 0 \rangle = 0 \nonumber$ . Since |0 should have the minimum energy, the L0-decomposition of V is then $\begin{array} { r } { { \mathcal V } \ = \ \bigoplus _ { n = 0 } ^ { \infty } \mathcal V _ { n } } \end{array}$ where $\begin{array} { r l r } { \mathcal { V } _ { 0 } } & { { } = } & { \mathbb { C } | 0 \rangle } \end{array}$ . Each space $\gamma _ { n }$ turns out to be finite dimensional, and we can think of $L _ { 0 }$ as defining a Z -grading on state space  .

The most important vertex operator in the theory is the stress-energy tensor T (z). The corresponding state is called the conformal vector ω: $Y ( \omega , z ) = T ( z )$ . This means that ω has modes $\omega _ { n } = L _ { n - 1 }$ that form a representation (4) of the Virasoro algebra . (This is Virthe algebraic expression for the requirement of conformal symmetry.) The conformal vector has energy 2: $\omega \in \mathcal { V } _ { 2 }$ .

So far our theory is seriously underdetermined. The most important axiom to help us to pin it down further is locality. With a little work, one can show that this reduces to the condition that the commutator $[ Y ( u , z ) , Y ( \nu , w ) ]$ of two vertex operators should be a finite linear combination of the Dirac delta $\delta ( z - w ) =$ $z ^ { - 1 } \Sigma _ { n = - \infty } ^ { \infty } ( w / z ) ^ { n }$ and its derivatives $( \partial ^ { k } / \partial w ^ { k } ) \delta ( z$ − w). Now, $( z - w ) ^ { k + 1 } ( \partial ^ { k } / \partial w ^ { k } ) \delta ( z - w ) \ = \ 0$ . To see this, look at the case $k = 1 \colon$ :

$$
\begin{array}{l} (z - w) ^ {2} \frac {\partial}{\partial w} \delta (z - w) \\ = \sum_ {n = - \infty} ^ {\infty} \left(n w ^ {n - 1} z ^ {- n + 1} - 2 n w ^ {n} z ^ {- n} + n w ^ {n + 1} z ^ {- n - 1}\right) \\ = \sum_ {n = - \infty} ^ {\infty} ((n + 1) - 2 n + (n - 1)) w ^ {n} z ^ {- n} = 0. \\ \end{array}
$$

The proof for general k is similar. Therefore, locality can be recast in an equivalent form as follows: given any $u , \upsilon \in \mathcal { V }$ , there is a positive number N such that

$$
(z - w) ^ {N} [ Y (u, z), Y (v, w) ] = 0. \tag {7}
$$

This equation may look strange. Why can we not simply divide out the $( z \mathrm { ~ - ~ } n ) ^ { N }$ and get that all vertex operators commute? The reason is that when formal power series are involved, there can be zero divisors. For example, it is easy to check that $\begin{array} { r l } { ( z - 1 ) \sum _ { n \in \mathbb { Z } } z ^ { n } = } & { { } } \end{array}$ 0. Locality in the form (7) is at the heart of VOAs; for instance, one can express it as a triply infinite sequence of identities that the modes must obey, and this emphasizes just how restrictive a condition it ${ \mathrm { i } } s ,$ and how correspondingly interesting it is to find examples of VOAs.

This completes the definition of a VOA. A consequence of these properties is that the modes $u _ { n }$ respect the L0-grading that we mentioned earlier. This means that if u has energy k and v has energy l, then $u _ { n } ( \nu )$ has energy $k + l - n - 1$ . The definition followed here is sometimes called a VOA of CFT-type, for obvious reasons. Sometimes in the literature some of these conditions are weakened or dropped. For example, much of the theory is independent of the existence of the conformal vector ω, although to us it will be crucial, for reasons that will be explained in the next subsection.

A VOA is simultaneously a physical and a mathematical object. We have emphasized their physical origins in order to help explain the motivation for studying them. We know they should be valuable to mathematics, simply because CFT is, and indeed this is the case, as we shall see in section 4. But from a purely mathematical point of view, they might appear somewhat ad hoc, as though we had a list of mathematical ingredients and said to ourselves, “Let’s consider this, and then have some of these, oh, and perhaps one of those too, but with the following extra assumption: . . . .” Fortunately, there are more abstract formulations of VOAs that make them appear much less arbitrary as mathematical structures. For example, Huang has shown that they can be regarded as “two-dimensionalized” Lie algebras, in the following sense. If you want to keep track of the Lie brackets in an expression such as $[ a , [ [ b , c ] , d ] ]$ (which is important since the Lie bracket is not an associative operation), you can do so with the help of a binary tree, and in fact it is easy to formulate Lie algebras in the language of such trees. If one then replaces binary trees by diagrams made out of spheres with legs, as we did with Feynman diagrams earlier, one obtains a structure that is equivalent to a VOA. (Of course, this is very far from a full explanation of what Huang did: his proof is extremely long.)

# 3.2 Basic Properties

We see from the definition sketched in the last subsection that a VOA is an infinite-dimensional Z -graded vector space with infinitely many products (namely u $* _ { n } \nu \ = \ u _ { n } ( \nu ) )$ , which obey infinitely many identities. Needless to say, it is not an easy definition, and there are no easy examples.

However, if we ignore the conformal symmetry (i.e., the conformal vector ω), then there are some simple, though uninteresting, examples. The easiest is the one-dimensional algebra ${ \mathcal { V } \mathrm { ~  ~ { ~ \mathbb ~ { ~ C ~ } ~ } ~ } } \mathbb { O } \rangle$ . More generally, a VOA  that obeys (7) with $N ~ = ~ 0$ is a commutative associative algebra with a unit $\mathrm { ~ 1 ~ } = \mathrm { ~ } | 0 \rangle$ . It also has a derivation $T \ = \ L _ { - 1 }$ , with respect to the product u $* \nu ~ = ~ u _ { - 1 } ( \nu ) \ d t$ : this means a linear map that obeys the product rule satisfied by derivatives, namely $T ( u * \nu ) = ( T u ) * \nu + u * ( T \nu )$ . The converse of this statement is true too: any such algebra is a VOA that obeys (7) with $N = 0 .$ . In these simple examples, the role of the derivation T is to recover the z-dependence of the vertex operator.

Therefore, we need N not to be zero in (7) if we want interesting examples. Likewise, the vertex operators $Y ( u , z )$ must be distributions (that is, they must involve doubly infinite sums) or again the VOA reduces to a commutative associative algebra.

It is also easy to show that in any VOA (again the existence of the conformal vector is not needed), the space $\mathcal { V } _ { 1 }$ is a Lie algebra, with Lie bracket given by $[ u \nu ] = u _ { 0 } ( \nu )$ . This is important because each $\gamma _ { n }$ will carry a representation of this Lie algebra, and $\mathcal { V } _ { 1 }$ generates continuous symmetries of the VOA (at least when $\mathcal { V } _ { 1 } ~ \neq ~ \{ 0 \} )$ . For a typical VOA V these Lie algebras are very familiar. For instance, for the VOAs associated with rational CFT, they are reductive, which means that they are a direct sum of copies of the trivial Lie algebra C with simple Lie algebras.

The existence of the conformal vector becomes important when one starts to consider the representation theory of VOAs. A -module is defined in a natural way. We shall not give full details here, but, roughly speaking, it is a space on which acts in such a way that as much as possible of the VOA structure is respected. For example,  will automatically be a module for itself, just as a group acts on itself in a simple way. (See representation theory [IV.9 §2] for an explanation of the latter.) A rational VOA is defined to be one that has the simplest representation theory: it has only finitely many irreducible -modules, and any V -module is a direct sum of irreducible ones. They are called rational VOAs because they are the VOAs that come from rational CFT. For these VOAs, acts irreducibly on itself.

Assume now that $\gamma$ is rational. Any irreducible $\mathcal { V } -$ module M will inherit from an $L _ { 0 }$ -grading by rational numbers, $M = \bigoplus _ { h } M _ { h }$ , into finite-dimensional spaces $M _ { h }$ . The character $\chi _ { M } ( \tau )$ is defined by

$$
\chi_ {M} (\tau) = \sum_ {h} \dim M _ {h} \mathrm{e} ^ {2 \pi \mathrm{i} \tau (h - c / 2 4)}, \tag {8}
$$

where c is the central charge. This definition arises naturally in CFT as well as in Lie theory (or affine Kac– Moody algebras), although the curious $" c / 2 4 , "$ needed for (9) below, is mysterious in Lie theory. (In CFT it has a natural explanation as a certain topological effect.) These characters converge for any τ in the upper halfplane H. They carry a representation of the modular group $\mathrm { S L } _ { 2 } ( \mathbb { Z } )$ :

$$
\chi_ {M} \left(\frac {a \tau + b}{c \tau + d}\right) = \sum_ {N \in \mathrm{Irr} (V)} \rho \left( \begin{array}{c c} a & b \\ c & d \end{array} \right) _ {M N} \chi_ {N} (\tau), \tag {9}
$$

where Irr(V ) denotes the (finite) set of irreducible V - modules, and $\rho ( { \bf \Pi } _ { c } ^ { a b } )$ is a matrix with complex entries, whose rows and columns are labeled by $M , N \in \operatorname { I r r } ( V )$ . Equation (9) holds for any $( \mathbf { \Pi } _ { c } ^ { a \ b } )$ in $\mathrm { S L } _ { 2 } ( Z )$ , i.e., for any integers $a , b , c ,$ d satisfying ad $b c = 1$ . The lengthy proof of (9), by Zhu, is perhaps the high point of VOA theory, and owes much to the intuitions of rational CFT. In the next section, we shall get some idea of why it is so important.

# 4 What Are VOAs Good For?

This section describes what are probably the two most significant applications of VOAs. But let us begin by listing (without any explanations) a few others. Inspired by the geometry of string theory, vertex operator (super)algebras have been assigned to manifolds, resulting in a powerful, though complicated, algebraic invariant of those manifolds that generalizes and enriches more classical data such as de Rham cohomology. VOAs associated with affine Kac–Moody algebras at “degenerate” levels k are deeply related to the geometric Langlands program. The modularity of affine algebra characters, as well as that of, for example, lattice theta functions, are all special cases of Zhu’s theorem, which places these modularities in a much broader context.

# 4.1 The Mathematical Formulation of CFT

Since the 1970s quantum field theory has had considerable success, especially in geometry, by studying classical structures using infinite-dimensional methods; this is a theme in particular of Atiyah’s school. Conformal field theories are a class of exceptionally symmetric quantum field theories, and they are also among the simplest nontrivial quantum field theories known. In the past two decades mathematics has feasted on this combination of symmetry and (relative) simplicity, often by “looping” or “complexifying” more classical structures, and the impact of CFT (or, equivalently, of string theory) has been especially significant and broad. In hindsight the importance of CFT to mathematics is not surprising: it is a coherent and intricate structure that straddles several disparate areas of mathematics, sprawling across geometry, number theory, analysis, combinatorics, and indeed algebra.

From this point of view, a crucial application of VOA theory has been to CFT itself. Quantum field theories are notoriously difficult to put on a rigorous mathematical footing. But the successful applications suggest that these difficulties are a symptom of mathematical profundity and subtlety rather than of irreparable mathematical incoherence. In this sense the situation is highly reminiscent of the deep conceptual challenges to eighteenth-century mathematicians that were raised by calculus. The definition of a VOA by Richard Borcherds makes the chiral algebra of a CFT completely rigorous, as well as concepts like the operator product expansion. Subsequent work (especially by Huang and Zhu) reconstructs from the VOA more and more of the CFT, in arbitrary genus. The resulting clarity makes the whole subject more accessible to, and hence exploitable by, mathematicians. Quantum field theories are here to stay in mathematics, and thanks to VOAs mathematicians are absorbing a large class of them completely and explicitly.

# 4.2 Monstrous Moonshine

In 1978 McKay noticed that 196 884 196 883. Why was this an interesting observation? Well, the number on the left is the first meaningful coefficient of the $j \cdot$ function [IV.1 §8]

$$
\begin{array}{l} j (\tau) = q ^ {- 1} + (7 4 4 +) 1 9 6 8 8 4 q + 2 1 4 9 3 7 6 0 q ^ {2} \\ + 8 6 4 2 9 9 9 7 0 q ^ {3} + \dots , \tag {10} \\ \end{array}
$$

the generator of all modular functions for ${ \mathrm { S L } } _ { 2 } ( \mathbb { Z } )$ . Recall that a modular function is a function $f ( \tau )$ that is meromorphic in the upper half-plane H and invariant under the usual action of $\mathrm { S L _ { 2 } ( Z ) }$ . It should also be meromorphic at the boundary points Q  i , which are called cusps; we did not mention this condition earlier. The j-function generates these functions in the sense that any such modular function $f ( \tau )$ can be written as a rational function poly $( j ( \tau ) ) / \mathrm { p o l y } ( j ( \tau ) )$ . In other words, j(τ) is a uniformizing function that identifies (H $\cup \mathbb { Q } \cup \{ { \mathrm { i } } { \infty } \} / { \mathrm { S L } } _ { 2 } ( \mathbb { Z } )$ with the Riemann sphere $\mathbb { C } \cup \infty$ . We bracketed the constant term 744 in (10) because although 744 was the traditional choice it can be freely replaced with any other number, including 0.

The number on the right in McKay’s observation is the dimension of the smallest nontrivial representation of the Monster, the most exceptional of the finite simple groups [V.7]. This relation between modular functions and the Monster was completely unexpected, as they seem to occupy completely independent spots in the mathematical universe. Conway, Norton, and others fleshed out and expanded McKay’s original observation by making a number of conjectures, collectively called Monstrous Moonshine. For instance, with every pair $( g , h )$ of commuting elements in the Monster (a group of size about $8 \times 1 0 ^ { 5 3 } )$ , we expect there to be associated a function $j _ { ( g , h ) } ( \tau )$ that generates all modular functions for some discrete subgroup $T _ { ( g , h ) }$ of $\mathrm { S L _ { 2 } ( Z ) }$ . The j-function would be assigned in the case $g = h = \mathrm { i d e n t i t y }$ .

The first major step toward proving these Moonshine conjectures was made by Frenkel, Lepowsky, and Meurman in the mid 1980s. They constructed an infinitedimensional vector space $V ^ { \natural }$ out of formal power series. They were motivated on the one hand by the vertex operators of string theory, and on the other by the formally similar distributions used in constructing affine algebra representations. This seemed a promising direction since for both string theory and affine algebra representations modular functions arise naturally. Together with a rich algebraic structure that came from these “vertex operators, $" V ^ { \natural }$ was also acted on in a natural way by the Monster group. Moreover, although $V ^ { \natural }$ is infinite dimensional, it comes packaged into finite-dimensional pieces $V ^ { \natural } = \bigoplus _ { n = - 1 } ^ { \infty } V _ { n } ^ { \natural }$ , and the “graded dimension” $) \textstyle \sum _ { n }$ dim $( V _ { n } ^ { \natural } ) q ^ { n }$ equals $j - 7 4 4$ . The action of the Monster sends each $V _ { n } ^ { \natural }$ to itself; that is, each space $V _ { n } ^ { \natural }$ itself carries a representation of the Monster. Frenkel, Lepowsky, and Meurman proposed that $V ^ { \natural }$ lies at the heart of the Monstrous Moonshine conjectures.

Borcherds was struck by the formal similarity between $V ^ { \natural }$ and the chiral algebras of CFTs, and by abstracting out their important algebraic properties he defined a new structure called a vertex (operator) algebra. His axioms clarified their relationship with (generalizations of) Kac–Moody algebras, and by 1992 he had proved the main Conway–Norton conjecture (which corresponds to the case where g is arbitrary but h is the identity in the conjecture given earlier). Although his definition of VOAs required a deep understanding of the physics of CFT, his elaborate proof of this Moonshine conjecture is purely algebraic.

We would now call $V ^ { \natural }$ a rational VOA with only one irreducible module (namely itself); its symmetry group is the Monster and its character (8) is $j ( \tau ) - 7 4 4$ . The removal of the constant term 744 from (10) is significant as it says that the Lie algebra $V _ { 1 } ^ { \natural }$ is trivial—this is necessary if the symmetry group is to be finite. It is conjectured that $V ^ { \natural }$ is the unique VOA with central charge $c \ = \ 2 4$ , trivial $\mathcal { V } _ { 1 }$ , and only one irreducible module. This is meant to be reminiscent of the leech lattice [I.4 §4], which is known to be the unique twenty-fourdimensional even self-dual lattice with no vectors of length $\sqrt { 2 }$ . Indeed, the Leech lattice plays a crucial role in the construction of $V ^ { \natural }$ .

Most of the Moonshine conjectures are still open and this deep connection between modular functions and the Monster is still somewhat mysterious. At the time of writing, however, VOAs still provide the only serious approach to the Moonshine conjectures.

Borcherds defined VOAs to clarify the chiral algebra of CFT and to tackle Monstrous Moonshine. For this work, he was awarded a Fields Medal in 1998.

# Further Reading

Borcherds, R. E. 1986. Vertex algebras, Kac–Moody algebras, and the Monster. Proceedings of the National Academy of Sciences of the USA 83:3068–71.   
. 1992. Monstrous Moonshine and monstrous Lie superalgebras. Inventiones Mathematicae 109:405–44.   
Di Francesco, P., P. Mathieu, and D. Sénéchal. 1996. Conformal Field Theory. New York: Springer.   
Gannon, T. 2006. Moonshine Beyond the Monster: The Bridge Connecting Algebra, Modular Forms and Physics. Cambridge: Cambridge University Press.   
Kac, V. G. 1998. Vertex Algebras for Beginners, 2nd edn. Providence, RI: American Mathematical Society.   
Lepowsky, J., and H. Li. 2004. Introduction to Vertex Operator Algebras and their Representations. Boston, MA: Birkhäuser.

# IV.18 Enumerative and Algebraic Combinatorics

Doron Zeilberger

# 1 Introduction

Enumeration, otherwise known as counting, is the oldest mathematical subject, while algebraic combinatorics is one of the youngest. Some cynics claim that algebraic combinatorics is not really a new subject but just a new name given to enumerative combinatorics in order to enhance its (former) poor image, but algebraic combinatorics is in fact the synthesis of two opposing trends: abstraction of the concrete and concretization of the abstract. The former trend dominated the first half of the twentieth century, starting with Hilbert’s “theological” proof of the fundamental theorem of invariants, in which he showed by abstract means that certain invariants existed, but not how to find them. The latter trend is dominating contemporary mathematics, thanks to the omnipresence of The Mighty Computer.

The abstraction trend consists of the categorization, conceptualization, structuralization, and fancification (in short, “bourbakization” [VI.96]) of mathematics. Enumeration did not escape this trend, and in the hands of such giants as Gian-Carlo Rota and Richard Stanley in America and Marco Schützenberger and Dominique Foata in France, classical, enumerative combinatorics became more conceptual, structural, and algebraic. However, as algebraic combinatorics has established itself as a fully fledged and separate mathematical speciality, the more recent trend toward the explicit, concrete, and constructive has left its mark as well. It has revealed that many algebraic structures have hidden combinatorial underpinnings; the attempts to unearth these have led to many fascinating discoveries and unsolved problems.

# 1.1 Enumeration

The fundamental theorem of enumeration, independently discovered by several anonymous cave dwellers, states that

$$
| A | = \sum_ {a \in A} 1.
$$

In words: the number of elements in A is the sum over all elements of A of the constant function 1.

While this formula is still useful after all these years, enumerating specific finite sets is no longer considered mathematics. A genuine mathematical fact has to incorporate infinitely many facts, and the generic enumeration problem is to enumerate not just one set but all the sets in an infinite family.

To be precise, given an infinite sequence of sets $\{ A _ { n } \} _ { n = 0 } ^ { \infty }$ , where each set $A _ { n }$ consists of objects satisfying some combinatorial specifications that depend on the parameter $n ,$ , answer the question: How many elements does $A _ { n }$ have?

In a moment we shall look at some examples. But before we can learn how to answer this kind of question, let us consider a meta-question: What is an answer?

This was posed, and beautifully answered, by Herbert Wilf. To give some background to Wilf’s meta-answer, let us examine answers to some famous instances of enumeration questions.

In the list below, when we are given a set $A _ { n }$ (which will change from example to example), we shall write $a _ { n }$ instead of $\vert A _ { n } \vert$ . That is, $a _ { n }$ will stand for the number of elements of $A _ { n }$ .

(i) I Ching. If $A _ { n }$ is the set of all subsets of $\{ 1 , \ldots , n \}$ , then $a _ { n } = 2 ^ { n }$ .   
(ii) Rabbi Levi Ben Gerson. If $A _ { n }$ is the set of permutations [III.68] on $\{ 1 , \ldots , n \}$ , then $a _ { n } = n !$ .   
(iii) Catalan. If $A _ { n }$ is the set of legal bracketings with n opening brackets and n closing brackets, then $a _ { n } =$ $( 2 n ) ! / ( n + 1$ )!n!. (A legal bracketing is a sequence of n opening brackets and n closing brackets such that at no point in the sequence has the number of closing brackets exceeded the number of opening brackets. For instance, when n = 2 the legal bracketings are [ ][ ] and [ [ ] ].)   
(iv) leonardo of pisa [VI.6]. Let $A _ { n }$ be the set of finite sequences that consist only of 1s and 2s and that sum to n. (For example, when n 4 the possible sequences are 1111, 112, 121, 211, and 22.) In this case, we have three equivalent answers as follows.

(i)

$$
a _ {n} = \frac {1}{\sqrt {5}} \left(\left(\frac {1 + \sqrt {5}}{2}\right) ^ {n + 1} - \left(\frac {1 - \sqrt {5}}{2}\right) ^ {n + 1}\right).
$$

(ii)

$$
a _ {n} = \sum_ {k = 0} ^ {\lfloor n / 2 \rfloor} \binom {n - k} {k}.
$$

(iii) $a _ { n } = F _ { n + 1 }$ , where $F _ { n }$ is the sequence defined by the recurrence $F _ { n } = F _ { n - 1 } + F _ { n - 2 } ,$ , subject to the initial conditions $F _ { 0 } = 0 , F _ { 1 } = 1$ .

(v) cayley [VI.46]. If $A _ { n }$ is the set of labeled trees on n vertices, then $a _ { n } = n ^ { n - 2 } .$ . (A tree is a connected graph [III.34] without cycles, and it is labeled if the vertices have distinct names.)   
(vi) If $A _ { n }$ is the set of labeled simple graphs with n vertices, then $a _ { n } = 2 ^ { n ( n - 1 ) / 2 }$ . (A graph is simple if it has neither loops nor multiple edges.)   
(vii) If $A _ { n }$ is the set of labeled connected simple graphs on n vertices (that is, graphs for which every vertex can be reached from every other by a path), then $a _ { n }$ is n! times the coefficient of $x ^ { n }$ in the power-series expansion of

$$
\log \bigg (\sum_ {k = 0} ^ {\infty} \frac {2 ^ {k (k - 1) / 2}}{k !} x ^ {k} \bigg).
$$

(viii) If $A _ { n }$ is the set of Latin squares of size n $( n \times n$ matrices each of whose rows and columns is a permutation of $\{ 1 , \ldots , n \} )$ , then not even a good approximation for $a _ { n }$ is known.

In 1982, Wilf defined an answer as follows.

Definition. An answer is a polynomial-time algorithm (in n) for computing $a _ { n }$ .

Wilf arrived at this definition after he refereed a paper proposing a “formula” for the answer to question (viii), and realized that its “computational complexity” exceeds that of the caveman’s formula of direct counting.

What is a “formula”? It is really an algorithm that inputs n and outputs $\scriptstyle { a _ { n } } .$ . For example, $a _ { n } ~ = ~ 2 ^ { n }$ is shorthand for the recursive algorithm

$$
i f n = 0 \text { then } a _ {n} = 1,
$$

$$
e l s e a _ {n} = 2 \cdot a _ {n - 1},
$$

which takes O(n) steps. However, using the algorithm

$$
i f n = 0 t h e n a _ {n} = 1,
$$

$$
e l s e \quad i f n \text {   is   odd,   then   } a _ {n} = 2 a _ {n - 1},
$$

$$
e l s e a _ {n} = a _ {n / 2} ^ {2}
$$

takes O(log n) steps, much faster than Wilf demands. In other cases, like enumerating self-avoiding walks, the best algorithm known is exponential, $O ( c ^ { n } )$ , and any lowering of the constant c is a major advance. (A self-avoiding walk is a sequence of points x0, $\pmb { x } _ { 1 } , \ldots , \pmb { x } _ { n }$ in the two-dimensional integer lattice, where each $\pmb { x } _ { i }$ is one of the four neighbors of $\pmb { x } _ { i - 1 }$ and no two of the $\pmb { x } _ { i }$ are equal.) Notwithstanding these exceptions, Wilf’s meta-answer is a very useful general guideline for evaluating answers.

Traditionally, the main customers of enumeration were probability and statistics. In fact, discrete probability is almost synonymous with enumerative combinatorics, since the probability of an event E occurring is the ratio of the number of successful cases divided by the total number. Also, statistical physics is, by and large, weighted enumeration of lattice models (see phase transitions and universality [IV.25]). About fifty years ago, another important customer came along: computer science. Here one is interested in the computational complexity [IV.20] of algorithms: that is, in the number of steps it takes to execute them.

# 2 Methods

The following tools are indispensable to the enumerative combinatorialist.

# 2.1 Decomposition

$$
| A \cup B | = | A | + | B | \quad (\text { if } A \cap B = \varnothing).
$$

In words: the size of the union of two disjoint sets equals the sum of their sizes.

$$
\left| A \times B \right| = \left| A \right| \cdot \left| B \right|.
$$

In words: the size of the Cartesian product of two sets (that is, the set of all pairs (a, b), where $a \in A$ and $b \in B )$ equals the product of their sizes.

$$
\left| A ^ {B} \right| = \left| A \right| ^ {| B |}.
$$

In words: the size of the set of functions from B to A equals the size of A raised to the power the size of B. For example, the number of 0–1 sequences of length n, which can be viewed as functions from $\{ 1 , 2 , \ldots , n \}$ to $\{ 0 , 1 \}$ , equals 2n.

# 2.2 Refinement

If

$$
A _ {n} = \bigcup_ {k} B _ {n k} \quad (\text { disjoint   union }),
$$

and if $b _ { n k } .$ , the number of elements of $B _ { n k } .$ , is “nice” (and even if it is not), then

$$
a _ {n} = \sum_ {k} b _ {n k}.
$$

The idea here is that it may be possible to take a set $A _ { n }$ that is difficult to count, and split it up into disjoint sets $B _ { n k }$ that are easier to count. For example, consider the set $A _ { n }$ of example (iv). This can be split into a disjoint union of subsets $B _ { n k }$ , where each $B _ { n k }$ consists of the sequences in $A _ { n }$ that have exactly k 2s. If there are k 2s, then there must be $n - 2 k \ 1 \mathrm { s } .$ , so $b _ { n k } = { \binom { n - k } { k } }$  n−k  . This yields answer (ii).

# 2.3 Recursion

Suppose that $A _ { n }$ can be decomposed in such a way that it is a combination of fundamental operations applied to the sets $A _ { n - 1 } , A _ { n - 2 } , \ldots , A _ { 0 }$ . Then $a _ { n }$ satisfies a recurrence relation of the form

$$
a _ {n} = P (a _ {n - 1}, a _ {n - 2}, \dots , a _ {0}).
$$

For example, let $A _ { n }$ be the set of example (iv). If a sequence in $A _ { n }$ starts with a 1, then the rest of the sequence must add up to $n - 1$ , and if it starts with a 2, then the rest must add up to n 2. Since when n $\geqslant 2$ exactly one of these possibilities occurs and both are possible, we can decompose $A _ { n }$ into $1 A _ { n - 1 }$ and $2 A _ { n - }$ −2, where $1 A _ { n - 1 }$ is shorthand for the set of all sequences that begin with a 1 and continue with a sequence in $A _ { n - 1 } ,$ , and $2 A _ { n - 2 }$ is defined similarly. Since the sizes of $1 A _ { n - 1 }$ and $2 A _ { n - 2 }$ are clearly $\scriptstyle a _ { n - 1 }$ and $\scriptstyle a _ { n - 2 } $ , it follows that $a _ { n } = a _ { n - 1 } + a _ { n - 2 }$ , which yields answer (iii).

If $A _ { n }$ is the set of legal bracketings with n pairs (example (iii)), then a typical legal bracketing can be written recursively as $[ L _ { 1 } ] L _ { 2 }$ , where $L _ { 1 }$ and $L _ { 2 }$ are smaller (possibly empty) legal bracketings. For example, if the bracketing is [ [ ] [ ] ] [ [ ] ] [ [ ] [ [ ] ] ] then $L _ { 1 } = \left[ \begin{array} { l l l } \end{array} \right]$ and $L _ { 2 } = [ [ ] ] [ [ ] [ [ ] ] ] . \operatorname { I f } L _ { 1 }$ has k pairs, then $L _ { 2 }$ has $n - 1 - k$ pairs. It follows that $A _ { n }$ can be identified with the union $\cup _ { k = 0 } ^ { n - 1 } A _ { k } \times A _ { n - 1 - k } ,$ and, taking cardinalities, $\begin{array} { r } { a _ { n } = \sum _ { k = 0 } ^ { n - 1 } a _ { k } a _ { n - 1 - k } } \end{array}$ − −. This is a nonlinear (in fact, quadratic) and nonlocal recurrence, but it is nevertheless one that satisfies Wilf’s dictum.

# 2.4 Generatingfunctionology

According to Wilf, who coined this neologism by making it the title of his classic book (a free download from his Web site, even though it is still in print!):

A generating function is a clothesline on which we hang up a sequence of numbers for display.

The method of generating functions is one of the most useful tools of the trade of enumeration. The generating function of a sequence, sometimes called its ztransform, is a discrete analogue of the laplace transform [III.91], and indeed goes back to laplace [VI.23] himself. If the sequence is $( a _ { n } ) _ { n = 0 } ^ { \infty } ,$ then its generating function $f ( x )$ is defined to be $\scriptstyle \sum _ { n = 0 } ^ { \infty } a _ { n } x ^ { n }$ . In other words, the terms of the sequence are regarded as the coefficients of a power series in x.

Generating functions are so useful because information about the sequence $( a _ { n } )$ translates to information about $f ( x )$ that is often easier to process, and after some manipulations one often gets additional information about $f ( x )$ that can be translated back into information about the sequence. For example, if $a _ { 0 } = a _ { 1 } = 1$ and $a _ { n } = a _ { n - 1 } + a _ { n - 2 }$ when $n \geqslant 2 ,$ , then we can do the following manipulations on f (x):

$$
\begin{array}{l} f (x) = \sum_ {n = 0} ^ {\infty} a _ {n} x ^ {n} = a _ {0} + a _ {1} x + \sum_ {n = 2} ^ {\infty} a _ {n} x ^ {n} \\ = 1 + x + \sum_ {n = 2} ^ {\infty} (a _ {n - 1} + a _ {n - 2}) x ^ {n} \\ = 1 + x + \sum_ {n = 2} ^ {\infty} a _ {n - 1} x ^ {n} + \sum_ {n = 2} ^ {\infty} a _ {n - 2} x ^ {n} \\ = 1 + x + x \sum_ {n = 2} ^ {\infty} a _ {n - 1} x ^ {n - 1} + x ^ {2} \sum_ {n = 2} ^ {\infty} a _ {n - 2} x ^ {n - 2} \\ = 1 + x + x (f (x) - 1) + x ^ {2} f (x) \\ = 1 + (x + x ^ {2}) f (x). \\ \end{array}
$$

It follows that

$$
f (x) = \frac {1}{1 - x - x ^ {2}}.
$$

If one performs a partial-fraction decomposition, and expands the two resulting terms in a Taylor series, then one can obtain answer (i) to example (iv).

# 3 Weight Enumeration

According to the modern approach, pioneered by Pólya, Tutte, and Schützenberger, generating functions are neither “generating,” nor are they functions. Rather, they are formal power series that are weight enumerators of combinatorial sets. (Usually, but not always, these sets are infinite: for a finite set the corresponding “power series” has only finitely many nonzero terms and is therefore a polynomial.)

A power series $\scriptstyle \sum _ { n = 0 } ^ { \infty } a _ { n } x ^ { n }$ is called formal when one sheds its analytical connotation as a Taylor series of a function, and thereby obviates the need to worry about convergence. For example, the sum $\scriptstyle \sum _ { n = 0 } n ! ^ { n ! } x ^ { n }$ is perfectly legal as a formal power series even though it converges only when $x = 0$ .

As for weight enumerators, consider the following situation. Suppose that we want to study the age distribution of a finite population. One way of doing this is to ask 121 questions. For each i between 0 and 120, we ask those whose age is i to raise their hand. Then we count each of these age-groups one by one, compiling a table of $a _ { i } ~ ( 0 \leqslant i \leqslant 1 2 0 )$ , and finally computing the generating function

$$
f (x) = \sum_ {i = 0} ^ {1 2 0} a _ {i} x ^ {i}.
$$

But if the size of the population is much less than 120, it is much more efficient, because fewer questions would be needed, to ask every person their age and then to declare the weight of a person of age i to be $x ^ { i }$ . Then the generating function is the sum of these weights. That is,

$$
f (x) = \sum_ {\text { persons }} x ^ {\text { age(person) }},
$$

which is a natural extension of the caveman’s formula of naive counting. Once we know $f ( x )$ we can easily compute statistically interesting quantities, like the average and the variance, which work out to be $\mu \ =$ $f ^ { \prime } ( 1 ) / f ( 1 )$ and $\sigma ^ { 2 } = f ^ { \prime \prime } ( 1 ) / f ( 1 ) + \mu - \mu ^ { 2 }$ , respectively.

The general scenario is that we have an interesting (finite or infinite) combinatorial set, let us call it A, and a certain numerical attribute, $\alpha : A  \mathbb { N } ,$ , which assigns to each element of A a natural number. (Here we allow 0 as a natural number.) Then the weight enumerator of A with respect to α is defined by the formula

$$
f (x) = \sum_ {a \in A} x ^ {\alpha (a)}.
$$

We shall also use the notation $| A | _ { x }$ for $f ( x )$ . Obviously, this equals

$$
\sum_ {n = 0} ^ {\infty} a _ {n} x ^ {n},
$$

where $a _ { n }$ is the number of members of A whose α equals n. Hence if we have some kind of explicit expression for $f ( x )$ , we immediately have an “explicit” expression for the actual sequence $a _ { n }$ assuming, that is, that one considers the operations needed to calculate the nth coefficient $a _ { n }$ of $f ( x )$ as constituting an explicit expression for $a _ { n }$ . Even if one does not, then it is still often possible to get a $\mathrm { ! \ " n i c e " }$ formula for $a _ { n }$ , or, failing this, to extract the asymptotics.

The fundamental operations for naive counting also hold for weighted counting: just replace | · | by | · |x. For example,

$$
\left| A \cup B \right| _ {x} = \left| A \right| _ {x} + \left| B \right| _ {x}
$$

$( { \mathrm { i f ~ } } A \cap B = \emptyset ) { \mathrm { ~ a n d } }$

$$
\left| A \times B \right| _ {x} = \left| A \right| _ {x} \cdot \left| B \right| _ {x}.
$$

Let us quickly see why the second of these is true. If the members of A and B are endowed with numerical attributes α and β, respectively, and one defines an attribute γ on A×B by letting $\gamma ( a , b )$ equal $\alpha ( a ) + \beta ( b )$ , then

$$
\begin{array}{l} | A \times B | _ {x} = \sum_ {(a, b) \in A \times B} x ^ {\gamma (a, b)} \\ = \sum_ {(a, b) \in A \times B} x ^ {\alpha (a) + \beta (b)} \\ = \sum_ {(a, b) \in A \times B} x ^ {\alpha (a)} \cdot x ^ {\beta (b)} \\ = \sum_ {a \in A} \sum_ {b \in B} x ^ {\alpha (a)} \cdot x ^ {\beta (b)} \\ = \left(\sum_ {a \in A} x ^ {\alpha (a)}\right) \cdot \left(\sum_ {b \in B} \cdot x ^ {\beta (b)}\right) \\ = | A | _ {x} \cdot | B | _ {x}. \\ \end{array}
$$

Let us see how these facts can be useful. First, consider the infinite set A, of all (finite) sequences of 1s and 2s, and let the attribute be “sum of entries.” Then the weight of 1221 is $x ^ { 6 } ,$ and, in general, the weight of a sequence $( a _ { 1 } \cdots a _ { r } )$ is $x ^ { a _ { 1 } + \cdots + a _ { k } }$ . The set A can be naturally decomposed as

$$
A = \{\phi \} \cup 1 A \cup 2 A,
$$

where $\phi$ is the empty word, and 1A is short for the set of all sequences obtained by prefixing a 1 to members of A, and analogously for 2A. Applying | · |x, we get

$$
| A | _ {x} = 1 + x | A | _ {x} + x ^ {2} | A | _ {x},
$$

which, in this simple case, can be solved explicitly, to yield, once again

$$
| A | _ {x} = \frac {1}{1 - x - x ^ {2}}.
$$

A legal bracketing L is either empty (in which case the weight is $x ^ { 0 } ~ = ~ 1 )$ , or else, as we have already noted, it can be written as $L = [ L _ { 1 } ] L _ { 2 }$ , where $L _ { 1 }$ and $L _ { 2 }$ are (shorter) legal bracketings. Conversely, whenever $L _ { 1 }$ and $L _ { 2 }$ are legal bracketings, so is $[ L _ { 1 } ] L _ { 2 }$ . Let $\mathcal { L }$ be the (infinite) set of all legal bracketings, and define the weight of a legal bracketing to be $x ^ { n }$ , where n is the number of bracket pairs [ ]. For example, the weight of [ ] is x and the weight of [ [ ] [ [ ] [ ] ] ] is $x ^ { 5 }$ . The set decomposes naturally as follows:

$$
\mathcal {L} = \{\phi \} \cup ([ \mathcal {L} ] \times \mathcal {L}),
$$

where φ denotes the empty word and $[ \mathcal { L } ] \times \mathcal { L }$ denotes the set of all words of the form $[ L _ { 1 } ] L _ { 2 }$ with $L _ { 1 }$ and $L _ { 2 }$ in L. This leads to the nonlinear (in fact, quadratic) equation

$$
| \mathcal {L} | _ {x} = 1 + x | \mathcal {L} | _ {x} ^ {2},
$$

which yields, thanks to the Babylonians, the explicit expression

$$
| \mathcal {L} | _ {x} = \frac {1 - \sqrt {1 - 4 x}}{2 x}.
$$

This in turn gives us the answer to example (iii) above, via Newton’s binomial theorem.

Legal bracketings are equivalent to so-called binary trees, that is, unlabeled ordered trees where every vertex has either no children or exactly two children. For instance, when we write the legal bracketing [ [ ] [ ] ] [ ] [ [ ] [ [ ] ] ] in the form $[ L _ { 1 } ] L _ { 2 }$ we can think of [ [ ] [ ] ] [ [ ] [ [ ] ] ] as the parent, with children $L _ { 1 } ~ = ~ [ ] [ ]$ and $L _ { 2 } ~ = ~ [ [ \bar { \bf \Phi } ] ] [ [ \bar { \bf \Phi } ] [ [ \bar { \bf \Phi } ] ] ]$ . Then $L _ { 1 } \ ' \mathbf { s }$ children are φ and [ ], while $L _ { 2 } '$ s are [ ] and [ [ ] [ [ ] ] ]. This process continues until we have reached φ down every branch of the family.

If we try to count penta-trees instead, where each vertex may only have exactly zero or five children, then the generating function, alias weight-enumerator, satisfies the quintic equation

$$
f = x + f ^ {5},
$$

which, according to abel [VI.33] and galois [VI.41], is not solvable by radicals (see the insolubility of the quintic [V.21]). However, solvability by radicals is not everything. More than 200 years ago, lagrange [VI.22] devised a beautiful and extremely useful formula for extracting the coefficients of the generating function from the equation it satisfies, now called the Lagrange inversion formula. Using it one can easily show that the number of complete k-ary trees with (k−1)m+1 leaves is

$$
\frac {(k m) !}{((k - 1) m + 1) ! m !}.
$$

A multivariate generalization of the Lagrange inversion formula, discovered by the great Bayesian probabilist I. J. Good, enables one to enumerate colored trees and many other extensions.

# 3.1 Enumeration Ansatzes

If one wants to turn enumerative combinatorics into a theory rather than a collection of solved problems, one needs to introduce classification, and enumeration paradigms for counting sequences. But since “paradigm” is such a pretentious word, let us use the much humbler German word “ansatz,” which roughly means “form of solution.”

Let $( a _ { n } ) _ { n = 0 } ^ { \infty }$ be a sequence, and let

$$
f (x) = \sum_ {n = 0} ^ {\infty} a _ {n} x ^ {n}
$$

be its generating function. If we know the “form” of $a _ { n }$ , we can often deduce the form of $f ( x )$ (and vice versa).

(i) If $a _ { n }$ is a polynomial in $n ,$ then $f ( x )$ has the form

$$
f (x) = \frac {P (x)}{(1 - x) ^ {d + 1}},
$$

where P is a polynomial function and d is the degree of the polynomial that describes $a _ { n }$ .

(ii) If $a _ { n }$ is a quasi-polynomial in n (i.e., there exists an integer N such that for each $r = 0 , \ldots , N - 1$ , the function $m \mapsto a _ { m N + r }$ is a polynomial in m), then, for some (finite) sequence of integers $d _ { 1 } , d _ { 2 } , \dots$ and some polynomial function $P ,$

$$
f (x) = \frac {P (x)}{(1 - x) ^ {d _ {1}} (1 - x ^ {2}) ^ {d _ {2}} (1 - x ^ {3}) ^ {d _ {3}} \dots}.
$$

(iii) If $a _ { n }$ is C-recursive, that is, if it satisfies a linear recurrence equation with constant coefficients

$$
a _ {n} = c _ {1} a _ {n - 1} + c _ {2} a _ {n - 2} + \dots + c _ {d} a _ {n - d}
$$

(a good example is the Fibonacci sequence), then $f ( x )$ is a rational function of $x :$ that is, $f ( x ) =$ $P ( x ) / Q ( x )$ , where P and $Q$ are polynomials.

(iv) If $a _ { n }$ satisfies a linear recurrence equation of the form

$$
\begin{array}{l} c _ {0} (n) a _ {n} = c _ {1} (n) a _ {n - 1} + c _ {2} (n) a _ {n - 2} \\ + \dots + c _ {d} (n) a _ {n - d}, \\ \end{array}
$$

where the coefficients $c _ { i } ( n )$ are polynomial in n, then it is said to be P-recursive. (For example, $a _ { n } =$ n! is P-recursive since we have the recurrence $a _ { n } =$ $n a _ { n - 1 } . )$ If this is the case, then $f ( x )$ is D-finite, which means that it satisfies a linear differential equation with polynomial coefficients (in x).

In the case of $a _ { n } =$ n! the recurrence $a _ { n } = n a _ { n - 1 }$ is first order. A natural example of a P-recursive sequence satisfying a higher-order linear recurrence with polynomial coefficients is the sequence that counts the number of involutions on $\{ 1 , \ldots , n \}$ . (An involution is a permutation that equals its inverse.) Let us call this number $w _ { n }$ . The sequence $( w _ { n } )$ satisfies the recurrence relation

$$
w _ {n} = w _ {n - 1} + (n - 1) w _ {n - 2}.
$$

This recurrence follows from the fact that in the permutation n belongs either to a 1-cycle or to a 2-cycle. The former case accounts for $w _ { n - 1 }$ of the involutions, and the latter for $( n - 1 ) w _ { n - 2 }$ of them. (There are $n - 1$ ways of choosing the cycle-mate, i, say, of n, and deleting the resulting cycle leaves an involution of the $n - 2$ elements $\{ 1 , \ldots , i - 1 , i + 1 , \ldots , n - 1 \} . )$

# 4 Bijective Methods

This last argument was a simple example of a bijective proof, in this case, of a recurrence for the number of involutions on n objects. Contrast it with the following proof.

The number of involutions of {1, . . . , n} with exactly k 2-cycles is

$$
\binom {n} {2 k} \frac {(2 k) !}{k ! 2 ^ {k}},
$$

because we must first choose the 2k elements that will participate in the k 2-cycles, and then match them up into (unordered) pairs, which can be done in

$$
(2 k - 1) (2 k - 3) \dots 1 = \frac {(2 k) !}{k ! 2 ^ {k}}
$$

ways. Hence

$$
w _ {n} = \sum_ {k} \binom{n}{2 k} \frac {(2 k) !}{k ! 2 ^ {k}}.
$$

Nowadays such sums can be handled completely automatically, and if one inputs this sum to the Maple package EKHAD (downloadable from my Web site), one would get the recurrence $w _ { n } \ = \ w _ { n - 1 } + ( n - 1 ) w _ { n - 2 }$ as the output, together with a (completely rigorous!) proof. While the so-called Wilf–Zeilberger (WZ) method is able to handle many such problems, there are many other cases where one still needs a human proof. In either case such proofs involve (algebraic, and sometimes analytic) manipulations. The great combinatorialist Adriano Garsia derogatorily calls such proofs “manipulatorics,” and real enumerators do not manipulate, or at least try to avoid it whenever possible. The preferred method of proof is by bijection [I.2 §2.2].

Suppose one has to prove that $\left| A _ { n } \right| = \left| B _ { n } \right|$ | for every $n ,$ where $A _ { n }$ and $B _ { n }$ are combinatorial families. The “ugly way” is to get, by some means or other, algebraic or analytic expressions for $a _ { n } = | A _ { n } |$ and $b _ { n } \ = \ | B _ { n } |$ . Then one manipulates $^ { a _ { n } , }$ getting another expression $a _ { n } ^ { \prime } ,$ which in turn leads to yet another expression $a _ { n } ^ { \prime \prime } ,$ and if one is patient enough, and clever enough, and in luck, or if the problem is not too deep, one eventually arrives at $b _ { n } ,$ , and the result follows.

On the other hand, the nice way of proving that $\left| A _ { n } \right| = \left| B _ { n } \right|$ is by constructing a (preferably nice) bijection $T _ { n } : A _ { n } \ \to \ B _ { n }$ , which immediately implies, as a corollary, that $\left| A _ { n } \right| = \left| B _ { n } \right|$ .

In addition to being more aesthetically pleasing, a bijective proof is also philosophically more satisfactory. In fact, the notion of (cardinal) number is a highly sophisticated derived notion based on the much more basic notion of being in bijection. Indeed, according to frege [VI.56], the cardinal numbers are equivalence classes, where the equivalence relation [I.2 §2.3] is “is in bijective correspondence with.” Saharon Shelah said that people have been exchanging objects, in a one-to-one way, since long before they started to count. Also, a bijective proof explains why the two sets are equinumerous, as opposed to just certifying the formal correctness of this fact.

For example, suppose that Noah had wanted to prove that there were as many male as female creatures in his Ark. One way of proving this would have been to count the males and count the females, and check that the two resulting numbers were indeed the same. But a much better, conceptual, proof would have been to note that there is an obvious one-to-one correspondence between the set M of males and the set F of females: the function $w : M \to F$ defined by w(x) = WifeOf(x) is a bijection, with inverse h : F → M defined by h(y) = HusbandOf(y).

A classic example of a bijective proof is Glaisher’s proof of euler’s [VI.19] “odd equals distinct” partition theorem. A partition of an integer n is a way of writing it as a sum of positive integers, where order does not matter. For example, 6 has eleven partitions: 6, 51, 42, 411, 33, 321, 3111, 222, 2211, 21111, 111111. (Here 3111 is shorthand for the sum 3 1 1 1, and so on. Since order does not matter, we count 3111 as the same partition of 6 as 1311, 1131, and 1113. It is convenient to write the partitions with their numbers in decreasing order, as we have done.)

A partition is called odd if all its parts are odd, and it is called distinct if all its parts are distinct. Let Odd(n) and Dis(n) be the sets of odd and distinct partitions of n, respectively. For example, Odd(6) = {51, 33, 3111, 111111} and Dis(6) = {6, 51, 42, 321}. Euler proved that |Odd(n)| = |Dis(n)| for all n. His “manipulatorics” proof goes as follows. Let o(n) and d(n) be the number of odd and distinct partitions of n, respectively, and let us define the generating functions

$$
f (q) = \sum_ {n = 0} ^ {\infty} o (n) q ^ {n} \quad \text { and } \quad g (q) = \sum_ {n = 0} ^ {\infty} d (n) q ^ {n}.
$$

With the help of the “multiplication principle” for weighted counting, Euler showed that

$$
f (q) = \prod_ {i = 0} ^ {\infty} \frac {1}{1 - q ^ {2 i + 1}} \quad \text { and } \quad g (q) = \prod_ {i = 0} ^ {\infty} (1 + q ^ {i}).
$$

Using the algebraic identity $1 + y = ( 1 - y ^ { 2 } ) / ( 1 - y )$ , we have

$$
\begin{array}{l} \prod_ {i = 0} ^ {\infty} (1 + q ^ {i}) = \prod_ {i = 0} ^ {\infty} \frac {1 - q ^ {2 i}}{1 - q ^ {i}} \\ = \frac {\prod_ {i = 0} ^ {\infty} (1 - q ^ {2 i})}{\prod_ {i = 0} ^ {\infty} (1 - q ^ {2 i}) \prod_ {i = 0} ^ {\infty} (1 - q ^ {2 i + 1})} \\ = \prod_ {i = 0} ^ {\infty} \frac {1}{1 - q ^ {2 i + 1}}. \\ \end{array}
$$

Hence $g ( q ) ~ = ~ f ( q )$ , and the identity $o ( n ) ~ = ~ d ( n )$ follows by extracting the coefficient of $q ^ { n }$ .

For a very long time, these kinds of manipulation were considered to belong to the realm of analysis, and in order to justify the manipulations of the infinite series and products, one talked about the “region of convergence,” usually |q| < 1, and every step had to be justified by the appropriate analytical theorem. Only relatively recently did people come to realize that no analysis need be involved: everything makes sense in the completely elementary and much more rigorous (from the philosophical viewpoint) algebra of formal power series. One still needs to worry about convergence, so as to exclude, for example, an infinite product like $\textstyle \prod _ { i = 0 } ^ { \infty } ( 1 + x )$ , but the notion of convergence in the ring of formal power series is much more user-friendly than its analytical namesake.

Even though invoking analysis was a red herring, Euler’s proof, while purely algebraic and elementary, is nevertheless still manipulatorics. It would be much nicer to find a direct bijection between the sets Dis(n) and Odd(n). Such a bijection was given by Glaisher. Given a distinct partition, write each of its parts as $2 ^ { r } \cdot s ,$ where s is odd, and replace it by $2 ^ { r }$ copies of s. (For example, $1 2 \ = \ 4 \cdot 3 ,$ so we would replace 12 by $3 + 3 + 3 + 3 . )$ The output is obviously a partition of the same integer n, but now into odd parts. For example, the partition (10, 5, 4) is transformed to the new partition (5, 5, 5, 1, 1, 1, 1). To define the inverse transformation, take an odd part a and count how many times it shows up. If it shows up m times, then write m in binary notation, $m = 2 ^ { s _ { 1 } } + \cdot \cdot \cdot + 2 ^ { s _ { k } }$ , and replace the m copies of a by the k parts: $2 ^ { s _ { 1 } } a , \ldots , 2 ^ { s _ { k } } a$ . It is not hard to check that if you do the first transformation to a partition in Dis(n) and then do the second transformation, you get back to the partition you started with.

When we perform algebraic (and logical, and even analytical) manipulations, we are really rearranging and combining symbols, and hence we are doing combinatorics in disguise. In fact, everything is combinatorics. All we need to do is to take the combinatorics out of the closet, and make it explicit. The plus sign turns into (disjoint) union, the multiplication sign becomes Cartesian product, and induction turns into recursion. But what about the combinatorial counterpart of the minus sign? In 1982, Garsia and Steven Milne filled this gap by producing an ingenious “involution principle” that enables one to translate the implication

$$
a = b \quad \text { and } \quad c = d \quad \Rightarrow \quad a - c = b - d
$$

into a bijective argument, in the sense that if $C \subset A$ and $D \subset B ,$ and there are natural bijections $f : A  B$ and $g : C \to D$ establishing that $| A | = | B |$ , and $| C | = | D |$ , then it is possible to construct an explicit bijection between $A \backslash C$ and $B \backslash D .$ . Let us define it in terms of people. Suppose that in a certain village all the adults are married, with the result that there is a natural bijection from the set of married men to the set of married women, m WifeOf(m), with its inverse $w \mapsto$ HusbandOf(w). In addition, some of the people have extramarital affairs, but only one per person, and all within the village. There is a natural bijection from the set of cheating men to the set of cheating women, called $m  \mathsf { M i s t r e s s O f } ( m )$ , with its inverse w  LoverOf(w). It follows that there are as many faithful men as there are faithful women. But how do we match them up? (One might imagine, for example, that each faithful man wants a faithful woman to go to church with him.)

Here is how it is done. A faithful man first asks his wife to come with him. If she is faithful, she agrees. If she is not, she has a lover, and that lover has a wife. So she tells her husband: “Sorry, hubby, I am going to the pub with my lover, but my lover’s wife may be free.” If this happens, then the man asks the wife of the lover of his wife to go with him, and if she is faithful, she agrees. If she is not he keeps asking the wife of the lover of the woman who has just rejected his proposal. Since the village is finite, he will eventually get to a faithful woman.

The reaction of the combinatorial enumeration community to the involution principle was mixed. On the one hand it had the universal appeal of a general principle, one that should be useful in many attempts to find bijective proofs of combinatorial identities. On the other hand, its universality is also a major drawback, since involution-principle proofs usually do not give any insight into the specific structures involved, and one feels a bit cheated. Such a proof answers the letter of the question, but it misses its spirit. Given a proof of this kind, one still hopes for a really natural, “involution-principle-free proof.” This is the case, for instance, with the celebrated Rogers–Ramanujan identity, which states that the number of partitions of an integer into parts that leave remainder 1 or 4 when divided by 5 equals the number of partitions of that integer with the property that the difference between any two parts is at least 2. For example, if $n = 7$ the cardinalities of 61, 4111, 1111111 and 7, 61, 52 are the same. Garsia and Milne invented their notorious principle in order to give a Rogers–Ramanujan bijection, thereby winning a \$50 prize from George Andrews. However, finding a really nice bijective proof is still an open problem.

A quintessential example of a bijective proof is Prüfer’s proof of cayley’s [VI.46] celebrated result that there are $n ^ { n - 2 }$ labeled trees on n vertices (example (v) earlier). Recall that a labeled tree is a labeled connected simple graph without cycles. Every tree has at least two vertices with only one neighbor (these are called leaves). A certain mapping called the Prüfer bijection associates with every labeled tree T a vector of integers $( a _ { 1 } , \ldots , a _ { n - 2 } )$ , with $1 \leqslant a _ { i } \leqslant n$ for each i. This vector is called its Prüfer code. Since there are $n ^ { n - 2 }$ such vectors, Cayley’s formula follows once we have defined the mapping $f$ : Trees  Codes and proved that it is indeed a bijection. This really needs four steps: defining $f ,$ defining its alleged inverse map $^ { g , }$ and proving that $g \circ f$ and $f \circ g$ are the identity maps on their respective domains.

The mapping $f$ is defined recursively as follows. If the tree has 2 vertices, then its code is the empty sequence. Otherwise, let $_ { a _ { 1 } }$ be the (sole) neighbor of the smallest leaf and let $( a _ { 2 } , \ldots , a _ { n - 2 } )$ be the code of the smaller tree obtained by deleting that leaf.

# 5 Exponential Generating Functions

So far, when we have discussed generating functions, we have been talking about ordinary generating functions (or OGFs). These are ideally suited for counting ordered structures like integer partitions, ordered trees, and words. But many combinatorial families are really sets, where the order is immaterial. For these the natural concept is that of an exponential generating function (or EGF).

The EGF of a sequence $\{ a ( n ) \} _ { n = 0 } ^ { \infty }$ is defined to be

$$
\sum_ {n = 0} ^ {\infty} \frac {a (n)}{n !} x ^ {n}.
$$

Labeled objects can be often viewed as sets of smaller irreducible objects. For example, a permutation is the disjoint union of cycles, a set partition is the disjoint union of nonempty sets, a (labeled) forest is the disjoint union of labeled trees, and so on.

Suppose that we have two combinatorial families A and B, and suppose that there are a(n) labeled objects of size n in the A family, and b(n) in the B family. We can construct a new set of labeled objects $C = A \times B ,$ where the labels are disjoint and distinct, and define the size of a pair to be the sum of the sizes of the components. We have

$$
c (n) = \sum_ {k = 0} ^ {n} \binom {n} {k} a (k) b (n - k),
$$

since we must

(i) decide the size of the first component, k (an integer between 0 and $^ { n ) , }$ , which forces the size of the second component to be $n - k$ ,   
(ii) decide which of the n labels go to the first component $( { \binom { n } { k } }$ ways), and   
(iii) pick the objects for each component from the A and B families, respectively, using the available labels $( a ( k ) b ( n - k )$ ways).

Multiplying both sides by $x ^ { n } /$ n! and summing from $n = 0$ to n = ∞ yields

$$
\begin{array}{l} \sum_ {n = 0} ^ {\infty} \frac {c (n)}{n !} x ^ {n} = \sum_ {n = 0} ^ {\infty} \sum_ {k = 0} ^ {n} \frac {a (k)}{k !} x ^ {k} \frac {b (n - k)}{(n - k) !} x ^ {n - k} \\ = \bigg (\sum_ {k = 0} ^ {\infty} \frac {a (k)}{k !} x ^ {k} \bigg) \bigg (\sum_ {n - k = 0} ^ {\infty} \frac {b (n - k)}{(n - k) !} x ^ {n - k} \bigg). \\ \end{array}
$$

Hence EGF(C) = EGF(A) EGF(B). Iterating, we get

$$
\operatorname{EGF} \left(A _ {1} \times A _ {2} \times \dots \times A _ {k}\right) = \operatorname{EGF} \left(A _ {1}\right) \dots \operatorname{EGF} \left(A _ {k}\right).
$$

In particular, if all the $A _ { i }$ are the same, we have that the EGF of ordered k-tuples, $A ^ { k } ,$ , equals $\operatorname { [ E G F } ( A ) ] ^ { k } .$ . But if “order does not matter,” then the EGF of k-sets of Aobjects is $[ \mathrm { E G F } ( A ) ] ^ { k } / k !$ , since there are exactly k! ways of arranging a k-set into an ordered array (since all labels are distinct, all these objects are different). Summing from $k = 0$ to k   we get the “fundamental theorem of exponential generating functions.”

If B is a labeled combinatorial family that can be viewed as sets of “connected components” that belong to a combinatorial family A, then

$$
\operatorname{EGF} (B) = \exp [ \operatorname{EGF} (A) ].
$$

This useful theorem was part of the physics folklore for many years, and was also implicit in many older combinatorial proofs. However, it was explicated only in the early 1970s. It was fully “categorized” by means of Joyal’s theory of species, which grew to be a beautiful theory of enumeration in the hands of the école Québecoise (the Labelle and Bergeron frères, Leroux, and others).

Here are some venerable examples. Let us try to find the EGF of set partitions. That is, let us try to figure out an expression for

$$
\sum_ {n = 0} ^ {\infty} \frac {b (n)}{n !} x ^ {n},
$$

where b(n) (so-called Bell numbers) denotes the number of set partitions of an n-element set.

Recall that a set partition of a set A is a set of pairwisedisjoint nonempty subsets of $A , \left\{ A _ { 1 } , \ldots , A _ { r } \right\}$ , such that the union of all the $A _ { i }$ equals A. For example, the set partitions of the 2-element set 1, 2 are 1 , 2 and 1, 2 .

The atomic objects in this example are nonempty sets. (We think of a set A as being the “trivial” partition of itself into just one set.) Let a(n) be the number of ways of partitioning a set of size n into one nonempty set. Clearly, when $n = 0$ this cannot be done, so $a ( 0 ) = 0$ . When n $\geqslant 1$ there is exactly one way of doing it, so the EGF of the sequence a(n) is

$$
A (x) = 0 + \sum_ {n = 1} ^ {\infty} \frac {1}{n !} x ^ {n} = \mathrm{e} ^ {x} - 1.
$$

It follows immediately from the fundamental theorem that

$$
\sum_ {n = 0} ^ {\infty} \frac {b (n)}{n !} x ^ {n} = \mathrm{e} ^ {\mathrm{e} ^ {x} - 1}, \tag {1}
$$

an identity of Bell. Nowadays, with computer algebra systems, this can be used immediately to crank out the first 100 terms of the sequence b(n). For example, in Maple one simply types

$$
\text { taylor } (\exp (\exp (x) - 1), x = 0, 1 0 1);
$$

so this is definitely an answer in the Wilfian sense. We can also easily derive recurrences (albeit ones that need at least O(n) memory), by differentiating both sides of (1) and comparing coefficients.

That was really easy, so let us go on and prove something much deeper. How about an EGF-style proof of Levi Ben Gerson’s celebrated formula for the number of permutations on n objects, n! (example (ii) earlier)? Every permutation can be decomposed into a disjoint union of cycles, so the atomic objects are now cycles. How many n-cycles are there? The answer is of course (n − 1)!, since $( a _ { 1 } , a _ { 2 } , \ldots , a _ { n } )$ is the same as $( a _ { 2 } , a _ { 3 } , \ldots , a _ { n } , a _ { 1 } )$ , which is the same as $( a _ { 3 } , \ldots , a _ { n } , a _ { 1 } , a _ { 2 } )$ , etc., which means that we can pick the first entry arbitrarily, after which we have (n  1)! choices for placing the remaining entries. The EGF for cycles is therefore

$$
\begin{array}{l} \sum_ {n = 1} ^ {\infty} \frac {(n - 1) !}{n !} x ^ {n} = \sum_ {n = 1} ^ {\infty} \frac {1}{n} x ^ {n} \\ = - \log (1 - x) = \log (1 - x) ^ {- 1}. \\ \end{array}
$$

Using the fundamental theorem, we get that the EGF of permutations is

$$
\exp (\log (1 - x) ^ {- 1}) = (1 - x) ^ {- 1} = \sum_ {n = 0} ^ {\infty} x ^ {n} = \sum_ {n = 0} ^ {\infty} \frac {n !}{n !} x ^ {n},
$$

and voilà we have a beautiful new proof that the number of permutations on n objects is n!.

This argument may not look very impressive. But a slight modification leads immediately to the (ordinary) generating function for the number of permutations on $\{ 1 , \ldots , n \}$ with exactly k cycles, which we shall denote by $c ( n , k )$ . Here we are fixing n and letting k vary, so the generating function is $\begin{array} { r } { C _ { n } ( \alpha ) = \sum _ { k = 0 } ^ { n } c ( n , k ) \alpha ^ { k } } \end{array}$ . All we have to do to calculate this is go from naive counting to weighted counting, and assign to each permutation the weight $\alpha ^ { \# \mathrm { c y c l e s } }$ . The fundamental theorem of exponential generating functions carries over word-forword to weighted counting. The weighted EGF for cycles is $\alpha \log ( 1 - x ) ^ { - 1 }$ , so the weighted EGF for permutations is

$$
\exp (\alpha \cdot \log (1 - x) ^ {- 1}) = (1 - x) ^ {- \alpha} = \sum_ {n = 0} ^ {\infty} \frac {(\alpha) _ {n}}{n !} x ^ {n},
$$

where

$$
(\alpha) _ {n} = \alpha (\alpha + 1) \cdot \cdot \cdot (\alpha + n - 1)
$$

is the so-called rising factorial. We have therefore derived the far less trivial result that the number of permutations of $\{ 1 , \ldots , n \}$ with exactly k cycles equals the coefficient of $\alpha ^ { k }$ in $( \alpha ) _ { n }$ .

About ten years ago (Ehrenpreis and Zeilberger 1994) I used this technique to give a combinatorial proof of the Pythagorean theorem in the form

$$
\sin^ {2} z + \cos^ {2} z = 1.
$$

The functions sin z and cos z are the weighted EGFs for increasing sequences of odd and even lengths, respectively, with weight $( - 1 ) ^ { [ \mathrm { l e n g t h } / 2 ] }$ . Hence the left-hand side is the weighted EGF for ordered pairs of increasing sequences

$$
a _ {1} <   \dots <   a _ {k}, \quad b _ {1} <   \dots <   b _ {r},
$$

such that k and r have the same parity, the sets $\{ a _ { 1 } , \ldots , a _ { k } \}$ and $\{ b _ { 1 } , \ldots , b _ { r } \}$ are disjoint, and the union of the two sets is $\{ 1 , 2 , \ldots , k + r \}$ . There is a killer involution on these sets of pairs defined as follows.

If $a _ { k } < b _ { r }$ then map the pair to

$$
a _ {1} <   \dots <   a _ {k} <   b _ {r}, \quad b _ {1} <   \dots <   b _ {r - 1}.
$$

and otherwise map it to

$$
a _ {1} <   \dots <   a _ {k - 1}, \quad b _ {1} <   \dots <   b _ {r} <   a _ {k}.
$$

For example, the pair

$$
\begin{array}{l l} 1, 3, 5, 6 & 2, 4, 7, 8, 9, 1 0, 1 1, 1 2, \end{array}
$$

whose sign is $( - 1 ) ^ { 2 } \cdot ( - 1 ) ^ { 4 } = 1$ , goes to the pair

$$
\begin{array}{l l} 1, 3, 5, 6, 1 2 & 2, 4, 7, 8, 9, 1 0, 1 1, \end{array}
$$

whose sign is $( - 1 ) ^ { 2 } \cdot ( - 1 ) ^ { 3 } = - 1$ (and vice versa).

Since this mapping changes the sign, and is an involution, all such pairs can be paired up into mutually canceling pairs. But this mapping is undefined for one special pair, namely the pair (empty, empty), whose weight is 1. Therefore, the EGF for the sum of the weights of all pairs is 1, which explains the right-hand side.

Yet another application of this method is a proof of André’s generating function for the number of up– down permutations. A permutation of $a _ { 1 } \cdots a _ { n }$ is called up–down (or sometimes zigzag) if $a _ { 1 } ~ < ~ a _ { 2 } ~ >$ $a _ { 3 } < a _ { 4 } > a _ { 5 } < \cdots$ . Let $a _ { n }$ be the number of up–down permutations. Then

$$
\sum_ {n = 0} ^ {\infty} \frac {a (n)}{n !} x ^ {n} = \sec x + \tan x.
$$

This is equivalent to saying that

$$
\cos x \cdot \left(\sum_ {n = 0} ^ {\infty} \frac {a (n)}{n !} x ^ {n}\right) = 1 + \sin x.
$$

Can you find the appropriate set and the killer involution?

# 6 Pólya–Redfield Enumeration

Often in enumeration it is easy enough to count labeled objects, but what about unlabeled ones? For example, the number of labeled (simple) graphs on n vertices (example (vi)) is trivially 2n(n−1)/2, but how many unlabeled graphs are there on n vertices? This is much harder, and in general there are no “nice” answers, but the best known way is via a powerful technique initiated by Pólya, which was largely anticipated by Redfield. Pólya enumeration lends itself very efficiently to counting chemical isomers, since, for example, all the carbon atoms “look the same.” Indeed, counting isomers was Pólya’s initial motivation (see mathematics and chemistry [VII.1 §2.3]).

The main idea is to view unlabeled objects as equivalence classes of easy-to-count labeled objects, and to count these equivalence classes. But what is the equivalence? The answer is that there is always a symmetry group [I.3 §2.1] involved, and it leads to a natural equivalence relation. Let the symmetry group be $G ,$ and let the set of labeled objects be A. Then two objects a and b of A are regarded as equivalent if $b = g ( a )$ for some member g of the group G. This means that there is some symmetry $_ g$ in the group G that transforms a to b. This is easily seen to be an equivalence relation and the equivalence classes are the sets

$$
\operatorname{Orbit} (a) = \{g (a) \mid g \in G \}, \quad a \in A,
$$

which are known as orbits. Calling each orbit a “family,” we have the task of counting the number of families. Note that G is a subgroup of the group of permutations of the finite set A.

Suppose that there is a picnic consisting of many families and we want to count the number of families. One way would be to define some “canonical head” of each family, say “mother,” and count the number of mothers. But some daughters look like mothers, so this is not so easy. On the other hand, you cannot just count everybody, since then you would count each family several times. The problem is that “naive” counting of people (or objects) is giving a credit of 1 to each person, and this is inappropriate if we are trying to count families. If instead we were to ask each person “How big is your family?” and add to our count the reciprocal of that number, then the calculation would come out just right, since a family of size k would get a credit of 1/k for each of its members, and would therefore have been counted exactly once by the end. Going back to counting orbits, we see by the same reasoning that their number is

$$
\sum_ {a \in A} \frac {1}{| \operatorname{Orbit} (a) |}.
$$

The conceptual opposite of “orbit of $\boldsymbol { a } ^ { \flat }$ is the subgroup of members of G that fix a:

$$
\operatorname{Fix} (a) = \{g \in G \mid g (a) = a \}.
$$

(This is sometimes known as the stabilizer of a.) To each element $b \ = \ g a$ in the orbit of a, we can associate the left coset g Fix(a) of Fix(a). This association turns out to be a well-defined one-to-one correspondence between the orbit of a and the cosets of Fix(a)

in G, from which it follows that the size of Orbit(a) is |G/ Fix(a)|. We can therefore substitute $\left| \operatorname { F i x } ( a ) \right| / \left| G \right|$ for 1/|Orbit(a)| in the previous formula, which implies that the number of orbits is

$$
\frac {1}{| G |} \sum_ {a \in A} | \operatorname{Fix} (a) |.
$$

Let us use the notation χ(statement) to stand for 1 if the statement is true and 0 if it is false. Then

$$
\begin{array}{l} \frac {1}{| G |} \sum_ {a \in A} | \operatorname{Fix} (a) | = \frac {1}{| G |} \sum_ {a \in A} \sum_ {g \in G} \chi (g (a) = a) \\ = \frac {1}{| G |} \sum_ {g \in G} \sum_ {a \in A} \chi (g (a) = a) \\ = \frac {1}{| G |} \sum_ {g \in G} \operatorname{fix} (g), \\ \end{array}
$$

where fix(g) is the number of fixed points of g (when g is viewed as a permutation of A). We have just proved what used to be called Burnside’s lemma, but it goes back to cauchy [VI.29] and frobenius [VI.58]. It states that the total number of orbits equals the average number of fixed points of $^ { g , }$ over all transformations g in G. If the group G is the full symmetric group of all the permutations of A, then the average number of fixed points equals 1 (since in this trivial case there is only one orbit!).

Enter Pólya. The objects that he was interested in counting (e.g., chemical isomers, or colorings of the faces of the cube) were all naturally functions from an underlying set to a set of colors (or atoms). Let us call the underlying set U and the set of colors C. A symmetry of U gives rise in a natural way to a transformation of the set of functions $f : U \to C$ . Given a function $f$ one defines a new function $g f$ by $g ( f ) ( u ) = f ( g ( u ) )$ . (If we think of f as a coloring, then gf is the new coloring that assigns to u the color that f assigned to g(u).) Now let us think about the number of fixed points of g in the set of C-colorings of U. Such a fixed point is a coloring f that equals gf : that is, $f ( u ) = f ( g u )$ for every u. But then $f ( u ) ~ = ~ f ( g u ) ~ = ~ f ( g ^ { 2 } u ) ~ = ~ \cdot ~ \cdot ~ .$ , which means that, given any cycle of $g , f$ must assign the same color to all members of that cycle. It follows that the number of fixed colorings of $^ g$ is $c ^ { \# \mathrm { c y c l e s } ( g ) }$ , where $c = | C |$ is the number of colors.

Applying Burnside’s lemma, we may deduce that the number of different colorings of U (up to G-equivalence) is

$$
\frac {1}{| G |} \sum_ {g \in G} c ^ {\# \text {cycles} (g)},
$$

since an equivalence class of colorings is simply an orbit of one of the colorings in that class.

Here is a simple application. How many necklaces (without a clasp) are there that consist of $p$ beads (where p is a prime) and that use a different colors? The underlying set is $\{ 0 , \ldots , p - 1 \}$ , and the symmetry group is $\mathbb { Z } _ { p }$ , the cyclic group of order $p .$ As usual, regard the elements of the symmetry group as permutations of the set of beads. Since p is a prime, there are p  1 elements of $\mathbb { Z } _ { p }$ with one cycle (of length p), and one element (the identity permutation) with p cycles (all of length 1). It follows that the number of necklaces is

$$
\frac {1}{p} ((p - 1) \cdot a + 1 \cdot a ^ {p}) = a + \frac {a ^ {p} - a}{p}.
$$

In particular, since this number is necessarily an integer, we get as a bonus a combinatorial proof of fermat’s little theorem [III.58]: that $\boldsymbol { a } ^ { p }$ a is always a multiple of $p .$ Perhaps one day there will be an equally nice combinatorial proof of Fermat’s last theorem. All one has to do is to prove that there is no bijection from the union of the set of straight necklaces of size n using x colors, and the set of such necklaces using y colors, to the set of necklaces using z colors (with $n > 2 ,$ , of course).

If one wants to keep track of how many beads there are of each color, one simply replaces straight counting by weighted counting, and $c ^ { \# \mathrm { c y c l e s } ( g ) }$ is replaced by

$$
(x _ {1} + \dots + x _ {c}) ^ {\alpha_ {1}} \cdot (x _ {1} ^ {2} + \dots + x _ {c} ^ {2}) ^ {\alpha_ {2}} \dots
$$

(assuming that $_ g$ has $\alpha _ { 1 }$ 1-cycles, $\alpha _ { 2 }$ 2-cycles, etc.). The resulting expression is the celebrated cycle-index polynomial.

# 6.1 The Principle of Inclusion–Exclusion and Möbius Inversion

Another pillar of enumeration is the principle of inclusion–exclusion (nicknamed PIE). Suppose that there are n sins, $s _ { 1 } , \ldots , s _ { n }$ , that a person may succumb to, and suppose that for each set of sins $S , A s$ is the set of people who have all the sins in S (and possibly others). Then the number of good people (without sins) is

$$
\sum_ {S} (- 1) ^ {| S |} | A _ {S} |.
$$

For example, if the set A is the set of all permutations π of $\{ 1 , \ldots , n \}$ and the ith sin is having $\pi [ i ] ~ = ~ i ,$ then $| A _ { S } | = ( n - | S | ) !$ , and we get that the number of derangements (permutations without fixed points) is

$$
\sum_ {k = 0} ^ {n} (- 1) ^ {k} \binom {n} {k} (n - k)! = n! \sum_ {k = 0} ^ {n} (- 1) ^ {k} \frac {1}{k !},
$$

which yields the answer: “closest integer to n!/e.” This is sometimes called the “umbrella problem”: if on a rainy day n absent-minded people go to a party and leave an umbrella by the door, and if on their departure they each take a random umbrella, then the probability that nobody ends up with the right umbrella is about 1/e.

The PIE is a special case of Möbius inversion on general partially ordered sets (posets) where the poset happens to be the Boolean lattice. This realization was published in a seminal paper by Rota (1964) and reprinted in his collected works. It is considered by many to be the big bang that started modern algebraic combinatorics. Möbius’s original inversion formula is recovered when the partially ordered set is N and the partial order is divisibility.

A contemporary account of enumeration from the “algebraic” point of view can be found in a marvelous two-volume set by Stanley (2000), which I strongly recommend.

# 7 Algebraic Combinatorics

So far I have described one of the routes to algebraic combinatorics: abstraction and conceptualization of classical enumeration. The other route, “concretization of the abstract,” is almost everywhere dense in mathematics, and cannot be described in a few pages. Let me quote from the preface of the excellent New Perspectives in Algebraic Combinatorics by Billera et al. (1999).

Algebraic combinatorics involves the use of techniques from algebra, topology, and geometry in the solution of combinatorial problems, or the use of combinatorial methods to attack problems in these areas. Problems amenable to the methods of algebraic combinatorics arise in these or other areas of mathematics or from diverse parts of applied mathematics. Because of this interplay with many fields of mathematics, algebraic combinatorics is an area in which a wide variety of ideas and methods come together.

# 7.1 Tableaux

An interesting class of objects that initially came up in group representation theory, but that turned out to be useful in many other areas—such as, for example, the theory of algorithms—are Young tableaux. They were first used by Reverend Alfred Young to construct explicit bases for the irreducible representations [IV.9 §2] of the symmetric group [III.68]. For any partition $\lambda = \lambda _ { 1 } \cdot \cdot \cdot \lambda _ { k }$ of n, a Young tableau of shape λ is an array of k left-justified rows with $\lambda _ { 1 }$ entries in the first row, λ2 entries in the second row, and so on, such that every row and every column is increasing, and the set of entries is $\{ 1 , 2 , \ldots , n \}$ . For example, there are two standard Young tableaux whose shape is 22,

<table><tr><td>1</td><td>2</td><td>1</td><td>3</td></tr><tr><td>3</td><td>4</td><td>2</td><td>4</td></tr></table>

and three of shape 31,

<table><tr><td>1</td><td>2</td><td>3</td><td></td><td>1</td><td>2</td><td>4</td><td></td><td>1</td><td>3</td><td>4</td><td></td></tr><tr><td>4</td><td></td><td></td><td></td><td>3</td><td></td><td></td><td></td><td>2</td><td></td><td></td><td>.</td></tr></table>

Let $f _ { \lambda }$ be the number of standard Young tableaux of shape λ. For example, for $n = 4 \colon f _ { 4 } = 1 , f _ { 3 1 } = 3 , f _ { 2 2 } =$ $2 , f _ { 2 1 1 } = 3 ,$ and $f _ { 1 1 1 1 } = 1$ . The sum of the squares of these numbers is $1 ^ { 2 } + 3 ^ { 2 } + 2 ^ { 2 } + 3 ^ { 2 } + 1 ^ { 2 } = 2 4 = 4 !$ .

The number $f _ { \lambda }$ is the dimension of the irreducible representation parametrized by λ. It follows by a result in representation theory [IV.9] known as Frobenius reciprocity that the same is true for all n. In other words,

$$
\sum_ {\lambda \vdash n} f _ {\lambda} ^ {2} = n!,
$$

a result known as the Young–Frobenius identity. A gorgeous bijective proof of this identity, which has many beautiful properties, was given by Gilbert Robinson and Craige Schensted and later extended by Donald Knuth, and is now known as the Robinson–Schensted– Knuth correspondence. It inputs a permutation π = $\pi _ { 1 } \pi _ { 2 } \cdots \pi _ { n }$ , and outputs a pair of Young tableaux of the same shape, thereby proving the identity.

Algebraic combinatorics is currently a very active field, and as mathematics is becoming more and more concrete, constructive, and algorithmic, there are going to be many more combinatorial structures discovered in all areas of mathematics (and science!) and this will guarantee that algebraic combinatorialists will stay very busy for a long time to come.

# Further Reading

Billera, L. J., A. Bjorner, C. Greene, R. E. Simion, and R P. Stanley, eds. 1999. New Perspectives in Algebraic Combinatorics. Cambridge: Cambridge University Press.   
Ehrenpreis, L., and D. Zeilberger. 1994. Two EZ proofs of sin2 z + cos2 z = 1. American Mathematical Monthly 101: 691.   
Rota, G.-C. 1964. On the foundations of combinatorial theory. I. Theory of Möbius functions. Zeitschrift für Wahrscheinlichkeitstheorie und Verwandte Gebiete 2:340– 68.   
Stanley, R. P. 2000. Enumerative Combinatorics, volumes 1 and 2. Cambridge: Cambridge University Press.

# IV.19 Extremal and Probabilistic Combinatorics

Noga Alon and Michael Krivelevich

# 1 Combinatorics: An Introduction

# 1.1 Examples

It is hard to give a rigorous definition of combinatorics. Instead, let us start with a few examples to illustrate what the area is about.

(i) In the course of an examination of friendship between children some fifty years ago, the Hungarian sociologist Sandor Szalai observed that among any group of about twenty children he checked he could always find four children any two of whom were friends, or else four children no two of whom were friends. Despite the temptation to try to draw sociological conclusions, Szalai realized that this might well be a mathematical phenomenon rather than a sociological one. Indeed, a brief discussion with the mathematicians Erd˝os, Turán, and Sós convinced him this was the case. If X is any set of size 18 or more, and R is some symmetric relation [I.2 §2.3] on X, then there is always a subset S of X of size 4 with the following property: either xRy for any two distinct elements $x , y$ of S, or xRy for no two distinct elements x, y of S. In this case, X is a set of children and R is the relation “is friends with.” This mathematical fact is a special case of Ramsey’s theorem, which was proved by the economist and mathematician Frank Plumpton Ramsey in 1930. Ramsey’s theorem led to the development of Ramsey theory, a branch of extremal combinatorics, which will be discussed in the next section.   
(ii) In 1916, Schur was studying fermat’s last theorem [V.10]. It is sometimes possible to prove that a Diophantine equation has no solutions by showing that it has no solutions mod p for some prime p. However, Schur proved that for every integer k and every sufficiently large prime p, there are three integers a, b, and $^ { c , }$ none of them congruent to 0 mod p, such that $a ^ { k } + b ^ { k }$ is congruent to $c ^ { k }$ . Although this is a result in number theory, it has a relatively simple and purely combinatorial proof, which is another example of the many applications of Ramsey theory.   
(iii) When studying the number of real zeros of random polynomials, littlewood [VI.79] and Offord investigated in 1943 the following problem. Let $z _ { 1 } , z _ { 2 } ,$ , $\cdots , z _ { n }$ be n not-necessarily-distinct complex numbers,

each of modulus at least 1. One can form $2 ^ { n }$ sums by taking some subset of these numbers and adding them together (with the convention that if one takes the empty set, then the sum is 0). Littlewood and Offord wanted to know how many of these sums there could conceivably be such that the difference between any two of them had modulus less than 1. When $n \ = \ 2$ the answer is easily seen to be at most 2. There are four sums: $0 , z _ { 1 } , z _ { 2 }$ , and $z _ { 1 } + z _ { 2 }$ . You cannot choose both of the first two or both of the last two or you will have a difference of $z _ { 1 }$ , which has modulus at least 1. Kleitman and Katona proved that in general the maximum is $\scriptstyle { \binom { n } { \lfloor n / 2 \rfloor } }$ . Notice that a simple construction proves that this maximum can be achieved. Indeed, let $z _ { 1 } = z _ { 2 } = \cdot \cdot \cdot = z _ { n }$ and choose all sums of precisely $\lfloor n / 2 \rfloor$ of them. There are $\scriptstyle { \binom { n } { \lfloor n / 2 \rfloor } }$ such sums and they are all equal. The proof that one cannot do better than this uses tools from another area of extremal combinatorics, where the basic objects studied are systems of finite sets.

(iv) Consider a school in which there are m teachers $T _ { 1 } , T _ { 2 } , \dots , T _ { m }$ and n classes $C _ { 1 } , C _ { 2 } , \ldots , C _ { n }$ . The teacher $T _ { i }$ has to teach the class $C _ { j }$ for a specified number $p _ { i j }$ of lessons. What is the minimum possible number of periods in a complete timetable? Let $d _ { i }$ denote the total number of lessons the teacher $T _ { i }$ has to teach, and let $c _ { j }$ denote the total number of lessons the class $C _ { j }$ has to be taught. Clearly, the number of periods required for a complete schedule is at least as big as any $d _ { i }$ or $c _ { j } ,$ and thus at least as big as the maximum of all these numbers, which we denote by d. It turns out that this obvious lower bound of d is also an upper bound: it is always possible to fit all the lessons that need to be taught into d periods. This is a consequence of König’s theorem, which is a basic result in graph theory. Suppose now that the situation is not so simple: for every teacher $T _ { i }$ and every class $C _ { j }$ there is some specified set of d periods in which the teaching has to take place. Can we always find a feasible timetable with these more complicated constraints? Recent breakthroughs from a subject known as list coloring of graphs imply that it is always possible.   
(v) Given a map with several countries represented, how many colors do you need if you want to color the countries without giving any two adjacent countries the same color? Here we assume that each country forms a connected region in the plane. Of course, at least four colors may be necessary: think of Belgium, France, Germany, and Luxembourg, out of which any two have a common border. The four-color theorem

[V.12], proved by Appel and Haken in 1976, asserts that you never need more than four colors. The study of this problem led to numerous interesting questions and results about graph coloring.

(vi) Let S be an arbitrary subset of the two-dimensional lattice $\mathbb { Z } ^ { 2 }$ . For any two finite subsets $A , B \subset \mathbb { Z }$ we can think of the Cartesian product $A \times B$ as a sort of “combinatorial rectangle.” This set has size $| A | \left| B \right|$ (where X denotes the size of a set X), and we can define an obvious notion of the density $d _ { S } ( A , B )$ of S in $A \times B$ by the formula $d _ { S } ( A , B ) = \vert S \cap ( A \times B ) \vert / \vert A \vert \vert B \vert$ , which measures what proportion of the elements of $A \times B$ belong to S. For each k, let $d ( S , k )$ be the largest possible value of $d _ { S } ( A , B )$ if $| A | = | B | = k$ . What can we say about $d ( S , k )$ as k tends to infinity? One might guess that almost any behavior is possible, but, remarkably, basic results in extremal graph theory (about the so-called Turán numbers of complete bipartite graphs) imply that d(S, k) must always tend to 0 or 1.

(vii) Suppose that n basketball teams compete in a tournament and any two teams play each other exactly once. The organizers wish to award k prizes at the end of the tournament. It would be embarrassing if there ended up being a team that had not won a prize despite beating all the teams that had won a prize. However, unlikely though it might sound, it is quite possible that this will be the case whatever k teams they choose, at least if n is large enough. To demonstrate this is easy if one uses the probabilistic method, which is one of the most powerful techniques in combinatorics. For any fixed $k ,$ and all sufficiently large n, if the results of all the games are chosen randomly (and uniformly and independently), then there is a very high probability that for any k teams there is another team that beats all of them. Probabilistic combinatorics, which is one of the most active areas in modern combinatorics, started with the realization that probabilistic reasoning often provides simple solutions to problems of this type, problems that are often very hard to solve in any other way.

(viii) If G is a finite group of n elements, and H is a subgroup of size k in G, then there are n/k left cosets and n/k right cosets of H. Is there always a set of $n / k$ elements of G that contains a single representative of each right coset and a single representative of each left coset? Hall’s theorem, a basic result in graph theory, implies that there is. In fact, if H	 is another subgroup of size k in G, then there is always a set of n/k elements of G that contains a single representative of each right coset of H and a single representative of each left coset of $H ^ { \prime }$ . This may sound like a result in group theory, but it is really a (simple) result in combinatorics.

# 1.2 Topics

The examples described above illustrate some of the main themes of combinatorics. The subject, sometimes also called discrete mathematics, is a branch of mathematics that focuses on the study of discrete objects (as opposed to continuous ones) and their properties. Although combinatorics is probably as old as the human ability to count, the field has experienced tremendous growth during the last fifty years and has matured into a thriving area with its own set of problems, approaches, and methodology.

The examples above suggest that combinatorics is a basic mathematical discipline that plays a crucial role in the development of many other mathematical areas. In this essay we discuss some of the main aspects of this modern field, focusing on extremal and probabilistic combinatorics. (An account of combinatorial problems with a rather different flavor can be found in algebraic and enumerative combinatorics [IV.18].) It is, of course, impossible to cover the area fully in such a short article. A detailed account of the subject can be found in Graham, Grötschel, and Lovász (1995). Our main intention is to give a glimpse of the topics, methods, and applications illustrated by representative examples. The topics we discuss include extremal graph theory, Ramsey theory, the extremal theory of set systems, combinatorial number theory, combinatorial geometry, random graphs, and probabilistic combinatorics. The methods applied in the area include combinatorial techniques, probabilistic methods, tools from linear algebra, spectral techniques, and topological methods. We also discuss the algorithmic aspects and some of the many fascinating open problems in the area.

# 2 Extremal Combinatorics

Extremal combinatorics deals with the problem of determining or estimating the maximum or minimum possible size of a collection of finite objects that satisfies certain requirements. Such problems are often related to other areas, including computer science, information theory, number theory, and geometry. This branch of combinatorics has developed spectacularly over the last few decades (see, for example, Bollobás (1978), Jukna (2001), and their many references).

# 2.1 Extremal Graph Theory

A graph [III.34] is one of the very basic combinatorial structures. It consists of a set of points, called vertices, some of which are linked by edges. One can represent a graph visually by drawing the vertices as points in the plane and the edges as lines (or curves). However, formally a graph is more abstract: it is just a set together with a collection of pairs taken from the set. More precisely, it consists of a set V , called the vertex set, and a set E, called the edge set; the elements of E (the edges) are sets of the form u, v , where u and v are distinct elements of V. If u, v is an edge, we say that u and v are adjacent. The degree d(v) of a vertex v is the number of vertices adjacent to it.

Here are a number of simple definitions associated with graphs that have emerged as important. A path of length k from u to v in G is a sequence of distinct vertices $u = \nu _ { 0 } , \nu _ { 1 } , \ldots , \nu _ { k } = \nu$ , where $\nu _ { i }$ and $\nu _ { i + 1 }$ are adjacent for all i < k. If $\upsilon _ { 0 } = \upsilon _ { k }$ (but all vertices $\nu _ { i }$ for $i < k$ are distinct), this is called a cycle of length k, and is usually denoted by $C _ { k } .$ A graph G is connected if for any two vertices u, v of G there is a path from u to v. A complete graph $K _ { r }$ is a graph with r vertices such that any two of them are adjacent. A subgraph of a graph G is a graph that contains some of the vertices of G and some of its edges. A clique in G is a set of vertices in G such that any two of them are adjacent. The maximum size of a clique in G is called the clique number of G. Similarly, an independent set in G is a set of vertices in G with no two of them adjacent, and the independence number of G is the maximum size of an independent set in it.

Extremal graph theory deals with quantitative connections between various parameters of a graph, such as its numbers of vertices and edges, its clique number, or its independence number. In many cases a certain optimization problem involving these parameters has to be solved (for example, determining how big one parameter can be if another one is at most some given size), and its optimal solutions are the extremal graphs for this problem. Many important optimization problems that do not explicitly mention graphs can be reformulated, using the definitions above, as problems about extremal graphs.

# 2.1.1 Graph Coloring

Let us return to the map-coloring example discussed in the introduction. To translate the problem into mathematics, we can describe the map-coloring problem in terms of a graph G, as follows. The vertices of G correspond to the countries on the map, and two vertices are connected by an edge in G if and only if the corresponding countries share a common border. It is not hard to show that one can draw such a graph in such a way that no two edges cross each other: such graphs are called planar. Conversely, any planar graph arises in this way. Therefore, our problem is equivalent to the following: if you want to color the vertices of a planar graph so that no two adjacent vertices receive the same color, then how many colors do you need? (One can make the problem yet more mathematical by removing the nonmathematical notion of color. For example, one can assign to each vertex a positive integer instead.) Such a coloring is called proper. In this language, the four-color theorem states that every planar graph can be properly colored with four colors.

Here is another example of a graph-coloring problem. Suppose we must schedule meetings of several parliament committees. We do not wish to have two committees meeting at the same time if some parliament member belongs to both, so how many sessions do we need?

Again we can model this situation by using a graph G. The vertices of G represent the committees, with two vertices adjacent if and only if the corresponding committees share a member. A schedule is a function f that assigns to each committee one of k time slots. More mathematically, we can think of it as just a function from V to the set 1, 2, . . . , k . Let us call a schedule valid if no two adjacent vertices are assigned the same number. This corresponds to no two committees being assigned the same time slot if they share a member. The question then becomes, “What is the minimal value of k for which a valid schedule exists?”

The answer is called the chromatic number of the graph G, denoted χ(G): it is the smallest number of colors in any proper coloring of G. Notice that a coloring of a graph G is proper if and only if for each color the set of vertices of that color is independent. Therefore, χ(G) can also be defined as the smallest number of independent sets into which it is possible to partition the vertices of G. A graph is called k-colorable if it admits a k-coloring, or, equivalently, if it can be partitioned into k independent sets. Thus, χ(G) is the minimum k for which G is k-colorable.

Two simple examples are in order. If G is a complete graph $K _ { n }$ on n vertices, then obviously in any coloring of G all vertices get distinct colors, and thus n colors are necessary. Of course, n colors are also sufficient, so $\chi ( K _ { n } ) ~ = ~ n$ . If G is a cycle $C _ { 2 n + 1 }$ on $2 n + 1$ vertices, then easy parity arguments show that at least three colors are needed, and three colors are enough: color the vertices along the cycle alternately by colors 1 and 2, and then color the last vertex by color 3. Thus, $\chi ( C _ { 2 n + 1 } ) = 3 .$ .

It is not hard to prove that G is 2-colorable if and only if it does not contain a cycle of odd length. Graphs that are 2-colorable are usually called bipartite, since they split into two parts, with all the edges going from one part to the other. The easy characterization ends here, and no simple criterion equivalent to k-colorability is available for $k \geqslant 3$ . This is related to the fact that for each fixed $k \geqslant 3$ the computational problem of deciding whether a given graph is k-colorable is NP-hard, a notion discussed in computational complexity [IV.20].

Coloring is one of the most fundamental notions of graph theory, as a huge array of problems in this field and in related areas like computer science and operations research can be formulated in terms of graph coloring. Finding an optimal coloring of a graph is known to be a very hard task, both theoretically and practically.

There are two simple yet fundamental lower bounds on the chromatic number. First, as every color class in a proper coloring of a graph G forms an independent set, it cannot be bigger than the independence number of $G ,$ which is denoted by $\alpha ( G )$ . Therefore, at least $\vert V ( G ) \vert / \alpha ( G )$ colors are necessary. Secondly, if G contains a clique of size $k ,$ then k colors are needed to color that clique alone, and thus $\chi ( G ) \geqslant k$ . This implies that $\chi ( G ) \geqslant \omega ( G )$ , where $\omega ( G )$ is the clique number of G.

What about upper bounds on the chromatic number? One of the simplest approaches to coloring a graph is to do it greedily: put the vertices in some order and color them one by one, assigning to each one the smallest positive integer that has not already been assigned to one of its neighbors. While the greedy algorithm can sometimes be very inefficient (for example, it can color bipartite graphs in an unbounded number of colors, even though two colors are sufficient), it often works quite well. Observe that when applying the greedy algorithm, a color given to a vertex v is at most one more than the number of the neighbors of v preceding it in the chosen order, and is thus at most $d ( \nu ) + 1$ , where d(v) is the degree of v in G. It follows that if $\Delta ( G )$ is the maximum degree of $G ,$ then the greedy algorithm uses at most $\Delta ( G ) + 1$ colors. Therefore $\chi ( G ) \leqslant \Delta ( G ) + 1$ . This bound is tight for complete graphs and odd cycles, and, as shown by Brooks in 1941, those are the only cases: if G is a graph of maximum degree $\Delta ,$ then $\chi ( G ) \leqslant \Delta$ unless G contains a clique $K _ { \Delta + 1 } ,$ or $\Delta \ : = \ : 2$ and G contains an odd cycle.

It is also possible to color the edges of a graph, rather than the vertices. In this case a proper coloring is defined to be one where no two edges that meet at a vertex are given the same color. The chromatic index of $G ,$ denoted by $\chi ^ { \prime } ( G )$ , is the minimum k for which G admits a proper edge-coloring with k colors. For example, if G is the complete graph $K _ { 2 n } ,$ , then $\chi ^ { \prime } ( G ) = 2 n - 1$ . This turns out to be equivalent to the fact that it is possible to organize a round-robin tournament with 2n teams and fit it into $2 n - 1$ rounds: just ask the manager of a soccer league. It is also not hard to show that $\chi ^ { \prime } ( K _ { 2 n - 1 } ) = 2 n - 1$ . Since in any proper edge-coloring of G all edges of G that are incident to a vertex v get distinct colors, the chromatic index is obviously at least as big as the maximum degree. Equality holds for bipartite graphs, as proved by König in 1931, which implies the existence of a complete timetable using d periods in the problem of teachers and classes discussed in the introduction.

Remarkably, this trivial lower bound of $\chi ^ { \prime } ( G ) \geqslant \Delta ( G )$ is very close to the true behavior of $\chi ^ { \prime } ( G ) . \mathrm { A }$ fundamental theorem of Vizing from 1964 states that $\chi ^ { \prime } ( G )$ is always equal either to the maximum degree $\Delta ( G )$ or to $\Delta ( G ) + 1$ . Thus, the chromatic index of G is much easier to approximate than its chromatic number.

# 2.1.2 Excluded Subgraphs

If a graph G has n vertices and contains no triangle (that is, three vertices all joined to each other) then how many edges can it contain? If n is even, then you can split the vertex set into two equal parts A and B of size $n / 2$ and join every vertex in A to every vertex in B. The resulting graph G contains no triangles and has $n ^ { 2 } / 4$ edges. Moreover, adding another edge will automatically create a triangle (in fact, several triangles). But is this the densest possible triangle-free graph? A hundred years ago the answer was shown to be yes by Mantel. (A similar theorem holds when n is odd, but now A and B must have nearly equal sizes $( n + 1 ) / 2$ and $( n - 1 ) / 2 . )$

Let us look at a more general problem, where the role of the triangle is played by an arbitrary graph. More precisely, let H be any graph, with m vertices, say, and when n  m let us define ex(n, H) to be the maximum possible number of edges in a graph with n vertices that does not contain H as a subgraph. (The notation “ex” stands for “exclude.”) The function ex(n, H) is usually called the Turán number of H, for reasons that will become clear, and finding good approximations for it has been a central problem in extremal graph theory.

What kind of examples of graphs that do not contain H can we think of? One observation that gets us started is that if H has chromatic number r , then it cannot be a subgraph of a graph G with chromatic number less than r . (Why not? Because a proper $( r - 1 )$ -coloring of G provides us with a proper (r  1)-coloring of any subgraph of G.) So a promising approach is to look for a graph G with n vertices, chromatic number $r - 1$ , and as many edges as possible. This is easy to find. Our constraint is that the vertices can be partitioned into r  1 independent sets. Once we have done that, we may as well include all edges between those sets. The result is a complete (r 1)-partite graph. A routine calculation shows that in order to maximize the number of edges, one should partition into sets that have sizes as nearly equal as possible. (For example, if n  10 and $r \ = 4 ,$ , then we would partition into three sets of sizes 3, 3, and 4.)

The graph that satisfies this condition is called the Turán graph $T _ { r - 1 } ( n )$ and the number of edges it contains is denoted by $t _ { r - 1 } ( n )$ . We have just argued that ex $( n , H ) \geqslant t _ { r - 1 } ( n )$ , which can be shown to be at least as big as $( 1 - 1 / ( r - 1 ) ) { \binom { n } { 2 } }$ .

Turán’s contribution to this area was to give an exact solution, in 1941, for the most important case, when H is the complete graph $K _ { r }$ on $r$ vertices. He proved that $\mathrm { e x } ( n , K _ { r } )$ is not just at least $t _ { r - 1 } ( n )$ , but is actually equal to $t _ { r - 1 } ( n )$ . Moreover, the only $K _ { r }$ -free graph with n vertices and $\mathrm { e x } ( n , K _ { r } )$ edges is the Turán graph $T _ { r - 1 } ( n )$ . Turán’s paper is generally considered the starting point of extremal graph theory.

Later, Erd˝os, Stone, and Simonovits extended Turán’s theorem by proving that the above simple lower bound for $\operatorname { e x } ( n , H )$ is asymptotically tight for any fixed H with chromatic number at least 3. That is, if r is the chromatic number of $H ,$ then the ratio of $\operatorname { e x } ( n , H )$ to $t _ { r - 1 } ( n )$ tends to 1 as n tends to infinity.

Thus, the function ex(n, H) is well-understood for all nonbipartite graphs. Bipartite graphs are rather different, because their Turán numbers are much smaller: if H is bipartite, then $\exp ( n , H ) / n ^ { 2 }$ tends to zero. Determining the asymptotics of ex(n, H) in this case remains a challenging open problem with many unsettled questions. Indeed, the full story is unknown even for the very simple case when H is a cycle. Partial results obtained so far use a variety of techniques from different fields, including probability theory, number theory, and algebraic geometry.

# 2.1.3 Matchings and Cycles

Let G be a graph. A matching in G is a collection of edges in G of which no two share a vertex. A matching M in G is called perfect if every vertex belongs to one of the edges in M. (The idea is that the edges determine $\mathrm { a \ ^ { * } m a t c h { \Sigma } ^ { * } }$ for each vertex: the match for x is the vertex y for which xy is an edge of M.) Of course, for G to have a perfect matching it must have an even number of vertices.

One of the best-known theorems in graph theory is Hall’s theorem, which provides a necessary and sufficient condition for the existence of a perfect matching in a bipartite graph. What kind of condition can this be? It is very easy to write down a trivial necessary condition, as follows. Let G be a bipartite graph with vertex sets A and B of equal size. (If they do not have equal size, then clearly there is no perfect matching.) Given any subset S of A, let N(S) denote the set of all vertices in B that are joined to at least one vertex in S. If there is to be a matching, then it must be possible to assign to each vertex in S a distinct “match,” so obviously N(S) must have at least as many elements as S. Hall’s theorem, proved in 1935, asserts that, remarkably, this obvious necessary condition is also sufficient. That is, if N(S) is always at least as big as S, then there will be a perfect matching. More generally, if A is smaller than B, then the same condition guarantees that one can find a matching that includes every vertex in A (but leaves some vertices in B unmatched).

There is a useful reformulation of Hall’s theorem in terms of set systems. Let $S _ { 1 } , S _ { 2 } , \ldots , S _ { n }$ be a collection of sets, and suppose that we would like to find a system of distinct representatives: that is, a sequence $x _ { 1 } , x _ { 2 } , \ldots , x _ { n }$ such that $x _ { i }$ is an element of $S _ { i }$ and no two of the $x _ { i }$ are the same. Obviously this cannot be done if the union of some k of the sets $S _ { i }$ has size less than k. Again, this obvious necessary condition is sufficient. It is not hard to show that this assertion is equivalent to Hall’s theorem: let S be the union of the $S _ { i }$ and define a bipartite graph with vertex sets $\{ 1 , 2 , \ldots , n \}$ and S, joining i to x if and only if $x \in S _ { i } .$ Then a matching that includes all of the set $\{ 1 , 2 , \ldots , n \}$ picks out a system of distinct representatives: $x _ { i }$ is the element of S that is matched with i.

Hall’s theorem can be applied to solve the problem of finding a system of representatives for the right and left cosets of a subgroup H, mentioned in section 1.1. Define a bipartite graph F, whose two sides (of size n/k each) are the left and right cosets of H. A left coset $g _ { 1 } H$ is connected by an edge of F to a right coset $H g _ { 2 }$ if they share a common element. It is not difficult to show that F satisfies the Hall condition, and hence it has a perfect matching M. Choosing for each edge $( g _ { i } H , H g _ { j } )$ of M a common element of $g _ { i } H$ and $H g _ { j }$ , we obtain the required family of representatives.

There is also a necessary and sufficient condition for the existence of a perfect matching in a general (not-necessarily-bipartite) graph G. This is a theorem of Tutte, which we shall not state here.

Recall that $C _ { k }$ denotes a cycle of length k. A cycle is a very basic graph structure, and, as one might expect, there are many extremal results concerning cycles.

Suppose that G is a connected graph with no cycles. If you pick a vertex and look at its neighbors and then the neighbors of its neighbors, and so on, you will see that it has a tree-like structure. Indeed, such graphs are called trees. An easy exercise shows that any tree with n vertices has exactly $n - 1$ edges. It follows that every graph G on n vertices with at least n edges has a cycle. If you want to guarantee that this cycle has certain extra properties, then you may need more edges. For example, the theorem of Mantel mentioned earlier implies that a graph G with n vertices and more than $n ^ { 2 } / 4$ edges contains a triangle $C _ { 3 } = K _ { 3 }$ . One can also prove that a graph $G \ : = \ : ( V , E )$ with $| E | > { \frac { 1 } { 2 } } k ( | V | - 1 )$ has a cycle of length longer than k (and this is in fact a sharp result).

A Hamilton cycle in a graph G is a cycle that visits every vertex of G. This term originated in a game, invented by hamilton [VI.37] in 1857, the objective of which was to complete a Hamilton cycle in the graph of the dodecahedron. A graph containing a Hamilton cycle is called Hamiltonian. This concept is strongly related to the well-known traveling salesman problem [VII.5 §2]: you are given a graph with positive weights assigned to the edges, and you must find a Hamilton cycle for which the sum of the weights of its edges is minimized. There are many sufficient criteria for a graph to be Hamiltonian, quite a few of which are based on the sequence of degrees. For example, Dirac proved in 1952 that a graph on $n \geqslant 3$ vertices all of whose degrees are at least $n / 2$ is Hamiltonian.

# 2.2 Ramsey Theory

Ramsey theory is a systematic study of the following general phenomenon. Surprisingly often, a large structure of a certain kind has to contain a fairly large highly organized substructure, even if the structure itself is completely arbitrary and apparently chaotic. As succinctly put by the mathematician T. S. Motzkin, “Complete disorder is impossible.” One might expect that the simple and very general form of this paradigm ensures that it has many diverse manifestations in different mathematical areas, and this is indeed the case. (One should, however, bear in mind that some natural statements of this kind are false for nonobvious reasons.)

A very simple statement, which can be regarded as a basic prototype for what follows, is the pigeonhole principle. This states that if a set X of n objects is colored with s colors, then there must be a subset of X of size at least $n / s$ that uses just one color. Such a subset is called monochromatic.

The situation becomes more interesting if the set X has some additional structure. It then becomes natural to ask for a monochromatic subset that keeps some of the structure of X. However, it also becomes much less obvious whether such a subset exists. Ramsey theory consists of problems and theorems of this general kind. Although several Ramsey-type theorems had appeared before, Ramsey theory is traditionally regarded as having started with Ramsey’s theorem, proved in 1930. Ramsey took as his set X the set of all the edges in a complete graph, and the monochromatic subset he obtained consisted of all the edges of some complete subgraph. A precise statement of his theorem is as follows. Let k and l be integers greater than 1. Then there exists an integer n such that, however you color the edges of the complete graph with n vertices, using the two colors red and blue, there will either be k vertices such that all edges between them are red or l vertices such that all edges between them are blue. That is, a sufficiently large complete graph colored with two colors contains a largish complete subgraph that is monochromatic. Let $R ( k , l )$ denote the minimum number n with this property. In this language, the observation of Szalai, mentioned in the introduction, is that $R ( 4 , 4 ) \leqslant 2 0$ (in fact, $R ( 4 , 4 ) = 1 8 )$ . Actually, Ramsey’s theorem was more general, in that he allowed any number of colors, and the objects colored could be r -tuples of elements rather than just pairs, as one has when coloring graphs. The exact computation of small Ramsey numbers turns out to be a notoriously difficult task: even the value of $R ( 5 , 5 )$ is unknown at present.

The second cornerstone of Ramsey theory was laid by Erd˝os and Szekeres, who in 1935 wrote a paper containing several important Ramsey-type results. In particular, they proved the recursion $R ( k , l ) \ \leqslant$ $R ( k - 1 , l ) + R ( k , l - 1 )$ . Combined with the easy boundary conditions $R ( 2 , l ) = l , R ( k , 2 ) = k ,$ , the recursion leads to the estimate $R ( k , l ) \leqslant { \binom { k + l - 2 } { k - 1 } }$ . In particular, for the so-called diagonal case $k = l$ we obtain $R ( k , k ) < 4 ^ { k }$ . Remarkably, no improvement in the exponent of the latter estimate has been found so far. That ${ \mathrm { i } } s ,$ nobody has found an upper bound of the form $C ^ { k }$ for some $C < 4$ . The best lower bound known, which we shall discuss in section 3.2, is roughly $R ( k , k ) \geqslant 2 ^ { k / 2 }$ , so there is a rather substantial gap.

Another Ramsey-type statement, proved by Erd˝os and Szekeres, is of a geometric nature. They showed that for every $n \geqslant 3$ there exists a positive integer N such that, given any configuration of N points in the plane in general position (i.e., no three of them are on a line), there are n that form a convex n-gon. (It is instructive to prove that if $n = 4$ then N can be taken to be 5.) There are several proofs of this theorem, some using the general Ramsey theorem. It is conjectured that the smallest value of N that will do in order to ensure a convex n-gon is $2 ^ { n - 2 } + 1$ .

The classic Erd˝os–Szekeres paper also contains the following Ramsey-type result: any sequence of $n ^ { 2 } + 1$ distinct numbers contains a monotone (increasing or decreasing) subsequence of length n 1.

This provides a quick lower bound of $\sqrt { n }$ for a well-known problem of Ulam, asking for the typical length of a longest increasing subsequence of a random sequence of length n. A detailed description of the distribution of this length has recently been given by Baik, Deift, and Johansson.

In 1927 van der Waerden proved what became known as van der Waerden’s theorem: for all positive integers k and r there exists an integer W such that for every coloring of the set of integers $\{ 1 , 2 , \ldots , W \}$ using r colors, one of the colors contains an arithmetic progression of length k. The minimum W for which this is true is denoted by $W ( k , r )$ . Van der Waerden’s bounds for $W ( k , r )$ are enormous: they grow like an Ackermanntype function. A new proof of his theorem was found by Shelah in 1987, and yet another proof was given by Gowers in 2000, while he was studying the (much deeper) “density version” of the theorem, which will be described in section 2.4. These recent proofs provided improved upper bounds for $W ( k , r )$ , but the bestknown lower bound for this number, which is only exponential in k for each fixed $r ,$ is much smaller.

Even before van der Waerden, Schur proved in 1916 that for any positive integer r there exists an integer S(r ) such that for every r -coloring of $\{ 1 , \ldots , S ( r ) \}$ one of the colors contains a solution of the equation $x + y = z .$ . The proof can be derived rather easily from the general Ramsey theorem. Schur applied this statement to prove the following result, mentioned in section 1.1: for every k and all sufficiently large primes $\nu ,$ the equation $a ^ { k } + b ^ { k } = c ^ { k }$ has a nontrivial solution in the integers modulo p. To prove this result, assume that $p \geqslant S ( k )$ and consider the field [I.3 §2.2] $\mathbb { Z } _ { p }$ of integers mod $p .$ The nonzero elements of $\mathbb { Z } _ { p }$ form a group [I.3 §2.1] under multiplication. Let H be the subgroup of this group consisting of all kth powers: that ${ \mathrm { i } } \mathbf { s } ,$ $H = \{ x ^ { k } : x \in \mathbb { Z } _ { p } ^ { * } \}$ . It is not hard to show that the index r of H is the highest common factor of k and p  1, and in particular is at most k. The partition of $\mathbb { Z } _ { p } ^ { \ast }$ into the cosets of H can be thought of as an r -coloring of Z∗p . By Schur’s theorem there exist $x , y , z \in \{ 1 , \ldots , p - 1 \}$ that all have the same color—that is, they all belong to the same coset of $H .$ In other words, there exists a residue $d \in \mathbb { Z } _ { p } ^ { * }$ such that $x \ : = \ : d a ^ { k } , \gamma \ : = \ : d b ^ { k } , z \ : = \ : d c ^ { k }$ , and $d a ^ { k } + { \dot { d } } b ^ { k } = d c ^ { k }$ modulo $p .$ . The desired result follows if we multiply both sides by $d ^ { - 1 }$ .

Many additional Ramsey-type results can be found in Graham, Rothschild, and Spencer (1990) or in Graham, Grötschel, and Lovász (1995, chapter 25).

# 2.3 Extremal Theory of Set Systems

Graphs are one of the fundamental structures studied by combinatorialists, but there are others too. An important branch of the subject is the study of set systems. Most often, these are simply collections of subsets of some n-element set. For example, the collection of all subsets of the set $\{ 1 , 2 , \ldots , n \}$ of size at most $n / 3$ is a good example of a set system. An extremal problem in this area is any problem where the aim is to determine, or estimate, the maximum number of sets there can be in a set system that satisfies certain conditions. For example, one of the first results in the area was proved by Sperner in 1928. He looked at the following question: how large a collection of subsets can one choose from an n-element set in such a way that no set from the collection is a subset of any other? A simple example of a set system satisfying this condition is the collection of all sets of size $r ,$ for some r . From this it immediately follows that we can obtain a collection as large as the largif n is even and mial coefficient, which is if n is odd. $\binom { n } { n / 2 }$ $\binom { n } { ( n + 1 ) / 2 }$

Sperner showed that this is indeed the maximum possible size of such a collection. This result supplies a quick solution to the real analogue of the problem of Littlewood and Offord described in section 1.1. Suppose that $x _ { 1 } , x _ { 2 } , \ldots , x _ { n }$ are n not-necessarily-distinct real numbers, each of modulus at least 1. A first observation is that we may assume that all the $x _ { i }$ are positive, since if we replace a negative $x _ { i }$ by $- x _ { i }$ (which is positive), then we end up with exactly the same set of sums, but shifted by $\cdot x _ { i } .$ (To see this, compare a sum that used to involve $x _ { i }$ with the corresponding sum that does not involve $- x _ { i } ,$ and vice versa.) But now, if A is a proper subset of $B ,$ then some $x _ { i }$ belongs to B and not to A, so

$$
\sum_ {i \in B} x _ {i} - \sum_ {i \in A} x _ {i} \geqslant x _ {i} \geqslant 1.
$$

Therefore, the total number of subset sums you can find with any two differing by less than 1 is at most $\scriptstyle { \binom { n } { \lfloor n / 2 \rfloor } }$ , by Sperner’s theorem.

A set system is called an intersecting family if any two sets in the system intersect. Since a set and its complement cannot both belong to an intersecting family of subsets of $\{ 1 , 2 , \ldots , n \}$ , we see immediately that such a family can have size at most $2 ^ { n - 1 }$ . Moreover, this bound is achieved by, for example, the collection of all sets that contain the element 1. But what happens if we fix a k and assume in addition that all our sets have size k? We may assume that $n \geqslant 2 k$ , as otherwise the solution is trivial. Erd˝os, Ko, and Rado proved that the maximum is  $\binom { n - 1 } { k - 1 }$ . Here is a beautiful proof discovered later by Katona. Suppose you arrange the elements randomly around a circle. Then there are n ways of choosing k elements that are consecutive in this arrangement, and it is quite easy to convince yourself that at most k of these can intersect (if $n \geqslant 2 k )$ . So out of these n sets of size $k ,$ only k of them can belong to any given intersecting family. Now it is also easy to show that every set has an equal chance of being one of these n sets, and this proves (by a simple double-counting argument) that the largest possible proportion of sets in the family is $k / n$ . Therefore, the family itself has size at most $( k / n ) { \binom { n } { k } }$ , kwhich equals n−1k−1 . The original proof of Erd˝os, Ko, and Rado is more complicated than this, but it is important because it introduced a technique known as compression, which was used to solve many other extremal problems.

Let n > 2k be two positive integers. Suppose that you wish to color all subsets of the set $\{ 1 , 2 , \ldots , n \}$ of size k in such a way that any two sets with the same color intersect each other. What is the smallest number of colors you can use? It is not difficult to see that $n { - } 2 k { + } 2$ colors suffice. Indeed, one color class can be the family of all subsets of $\{ 1 , 2 , \ldots , 2 k - 1 \}$ , which is clearly an intersecting family. And then, for each i such that $2 k \leqslant$ $i \ \leqslant \ n ,$ , you can take the family of all subsets whose largest element is i. There are $n - 2 k + 1$ such families, and any set of size k belongs either to one of them or to the first family. Therefore, $n - 2 k + 2$ colors are enough.

Kneser conjectured in 1955 that this bound was tight: in other words, that if you have fewer than n $2 k + 2$ colors then you will have to give the same color to some pair of disjoint sets. This conjecture was proved by Lovász in 1978. His proof is topological, and uses the Borsuk–Ulam theorem. Several simpler proofs have been found since, but they are all based on the topological idea in the first proof. Since Lovász’s breakthrough, topological arguments have become an important part of the armory of researchers in combinatorics.

# 2.4 Combinatorial Number Theory

Number theory is one of the oldest branches of mathematics. At its core are problems about integers, but a sophisticated array of techniques has been developed to deal with those problems, and these techniques have often themselves been the basis for further study (see, for example, algebraic numbers [IV.1], analytic number theory [IV.2], and arithmetic geometry [IV.5]). However, some problems in number theory have yielded to the methods of combinatorics. Some of these problems are extremal problems with a combinatorial flavor, while others are classical problems in number theory where the existence of a combinatorial solution has been quite surprising. We describe below a few examples. Many more can be found in chapter 20 of Graham, Grötschel, and Lovász (1995), in Nathanson (1996), and in Tao and Vu (2006).

A simple but important notion in the area is that of a sumset. If A and B are two sets of integers, or more generally are two subsets of an abelian group [I.3 §2.1], then the sumset $A + B$ is defined to be $\{ a + b : a \in$ A, $b \in B \}$ . For instance, if $A = \{ 1 , 3 \}$ and $B = \{ 5 , 6 , 1 2 \}$ , then $A + B = \{ 6 , 7 , 8 , 9 , 1 3 , 1 5 \}$ . There are many results relating the size and structure of $A + B$ to those of A and B. For example, the Cauchy–Davenport theorem, which has numerous applications in additive number theory, is the statement that if p is a prime, and $A ,$ B are two nonempty subsets of $\mathbb { Z } _ { p }$ , then the size of $A + B$ is at least the minimum of $p$ and $| A | + | B | - 1$ . (Equality occurs if A and B are arithmetic progressions with the same common difference.) cauchy [VI.29] proved this theorem in 1813, and applied it to give a new proof of a lemma that lagrange [VI.22] had proved as part of his well-known 1770 paper that shows that every positive integer is a sum of four squares. Davenport formulated the theorem as a discrete analogue of a related conjecture of Khinchin about densities of sums of sequences of integers. The proofs given by Cauchy and by Davenport are combinatorial, but there is also a more recent algebraic proof, based on some properties of roots of polynomials. The advantage of the latter is that it provides many variants that do not seem to follow from the combinatorial approach. For example, let us define A B to be the set of all a b such that a $\in A , b \in B ,$ and a b. Then the smallest possible size of A B, given the sizes of A and B, is the minimum of p and $| A | + | B | - 2$ . Further extensions can be found in Nathanson (1996) and in Tao and Vu (2006).

The theorem of van der Waerden mentioned in section 2.2 implies that, however you color the positive integers with some finite number r of colors, there must be some color that contains arithmetic progressions of every length. Erd˝os and Turán conjectured in 1936 that this always holds for the “most popular” color class. More precisely, they conjectured that for any positive integer k and for any real number $\epsilon > 0 ,$ , there is a positive integer $n _ { 0 }$ such that if $n \ > \ n _ { 0 } ,$ , any set of at least n positive integers between 1 and n contains a k-term arithmetic progression. (Setting $\epsilon = r ^ { - 1 }$ one can easily deduce van der Waerden’s theorem from this.) After several partial results, this conjecture was proved by Szemerédi in 1975. His deep proof is combinatorial, and applies techniques from Ramsey theory and extremal graph theory. Furstenberg gave another proof in 1977, based on techniques of ergodic theory [V.9]. In 2000 Gowers gave a new proof, combining combinatorial arguments with tools from analytic number theory. This proof supplied a much better quantitative estimate. A related very recent spectacular result of Green and Tao asserts that there are arbitrarily long arithmetic progressions of prime numbers. Their proof combines number-theoretic techniques with the ergodic theory approach. Erd˝os conjectured that any infinite sequence $n _ { i }$ for which the sum $\textstyle \sum _ { i } ( 1 / n _ { i } )$ diverges contains arbitrarily long arithmetic progressions. This conjecture would imply the theorem of Green and Tao.

# 2.5 Discrete Geometry

Let P be a set of points and let L be a set of lines in the plane. Let us define an incidence to be a pair $( p , \ell )$ , where p is a point in P ,  is a line in $L ,$ and the point p lies on the line . Suppose that P contains m distinct points and L contains n distinct lines. How many incidences can there be? This is a geometrical problem, but again it has a strong flavor of extremal combinatorics. As such, it is typical of the area known as discrete (or combinatorial) geometry.

Let us write I(m, n) for the maximum number of incidences there can be between m points and n lines. Szemerédi and Trotter determined the asymptotic behavior of this quantity, up to a constant factor, for all possible values of m and n. There are two absolute positive constants $c _ { 1 } , c _ { 2 }$ such that, for all m, n,

$$
\begin{array}{l} c _ {1} (m ^ {2 / 3} n ^ {2 / 3} + m + n) \leqslant I (m, n) \\ \leqslant c _ {2} \left(m ^ {2 / 3} n ^ {2 / 3} + m + n\right). \\ \end{array}
$$

If $m > n ^ { 2 } \ { \mathrm { o r } } \ n > m ^ { 2 }$ then one can establish the lower bound by taking all m points on a single line, or all n lines through a single point, respectively. In the harder cases when m and n are closer to each other, one can prove it by letting $P$ contain all the points of $\mathrm { ~ a ~ l ~ } \sqrt { m } \mathrm { ] }$ by $\lfloor \sqrt { m } \rfloor$ grid, and by taking the n most “popular” lines: that is, the n lines that contain the most points of P. Establishing the upper bound is more difficult. The most elegant proof of it is due to Székely, and is based on the fact that, however you draw a graph with m vertices and more than 4m edges, you must have many pairs of edges that cross each other. (This is a rather simple consequence of the famous Euler formula connecting the numbers of vertices, edges, and regions in any drawing of a planar graph.) To bound the number of incidences between a set of points P and a set of lines L in the plane, one considers the graph whose vertices are the points $P ,$ and whose edges are all segments between consecutive points along a line in L. The desired bound is obtained by observing that the number of crossings in this graph does not exceed the number of pairs of lines in L, and yet should be large if there are many incidences.

Similar ideas can be used to give a partial answer to the following question: if you take n points in the plane, how many pairs $( x , y )$ of these points can there be with the distance from x to y equal to 1? It is not surprising that the two problems are related: the number of such pairs is the number of incidences between the given n points and the n unit circles that are centered at these points. Here, however, there is a large gap between the best known upper bound, which is $c n ^ { 4 / 3 }$ for some absolute constant $^ { c , }$ and the best known lower bound, which is only $n ^ { 1 + c ^ { \prime } / }$ / log log n for some constant $c ^ { \prime } > 0 \ .$ .

A fundamental theorem of Helly asserts that if you have a finite family $\mathcal { F }$ of at least $d + 1$ convex sets in $\mathbb { R } ^ { d }$ , and if any d + 1 of them have a point in common, then all sets in the family have a common point. Now let us start with a weaker assumption: given any p of the sets, some $d + 1$ of those p sets have a point in common. (Here p is some integer greater than d  1.) Can one then find a set X of at most C points such that each set in  contains a point in X, with C a constant that depends on p but not on the number of convex sets in the family? This question was raised by Hadwiger and Debrunner in 1957 and solved by Kleitman and Alon in 1992. The proof combines a “fractional version” of Helly’s theorem with the duality of linear programming [III.84] and various additional geometric results. Unfortunately, it gives a very poor estimate for C: even in two dimensions and with $p = 4 \mathrm { i }$ t is not known what the best possible value of C is.

This is just a small sample of problems and results in discrete geometry. Such results have been applied extensively in computational geometry and in combinatorial optimization in recent decades. Two good books on the subject are Pach and Agarwal (1995) and Matoušek (2002).

# 2.6 Tools

Many of the basic results in extremal combinatorics were obtained mainly by ingenuity and detailed reasoning. However, the subject has grown out of this early stage: several deep tools have been developed that have been essential to much of the recent progress in the area. In this subsection, we include a very brief description of some of these tools.

Szemerédi’s regularity lemma is a result in graph theory that has numerous applications in various areas, including combinatorial number theory, computational complexity, and, mainly, extremal graph theory. The precise statement of the lemma, which can be found, for example, in Bollobás (1978), is somewhat technical. The rough statement is that the vertex set of any large graph can be partitioned into a constant number of pieces of nearly equal size, so that the bipartite graphs between most pairs of pieces behave like random bipartite graphs. The strength of this lemma is that it applies to any graph, providing a rough approximation of its structure that enables one to extract a lot of information about it. A typical application is that a graph with “few” triangles can be “well-approximated” by a graph with no triangles. More precisely, for any $\epsilon > 0$ there exists $\delta > 0$ such that if G is a graph with n vertices and at most $\delta n ^ { 3 }$ triangles, then one can remove at most $\epsilon n ^ { 2 }$ edges from G and make it triangle free. This innocentlooking statement turns out to imply the case $k = 3$ o f Szemerédi’s theorem that was mentioned earlier.

Tools from linear and multilinear algebra play an essential role in extremal combinatorics. The most fruitful technique of this kind, which is possibly also the simplest, is the so-called dimension argument. In its simplest form, the method can be described as follows. In order to bound the cardinality of a discrete structure A, one maps its elements to distinct vectors in a vector space [I.3 §2.3], and proves that those vectors are linearly independent. It then follows that the size of A is at most the dimension of the vector space in question. An early application of this argument was found by Larman, Rogers, and Seidel in 1977. They wanted to know how many points it was possible to find in $\mathbb { R } ^ { n }$ that determine at most two distinct differences. An example of such a system is the set of all points whose coordinates consist of $n - 2$ 0s and two 1s. Notice, however, that these points all lie in the hyperplane of points whose coordinates add up to 2. So this actually provides us with an example in $\mathbb { R } ^ { n - 1 }$ . Therefore, we have a simple lower bound of $n ( n + 1 ) / 2$ . Larman, Rogers, and Seidel matched this with an upper bound of $\left( n + 1 \right) \left( n + 4 \right) / 2$ . They did this by associating with each point of such a set a polynomial in n variables, and by showing that these polynomials are linearly independent and all lie in a space of dimension $\left( n + 1 \right) \left( n + 4 \right) / 2$ . This has been improved by Blokhuis to $\left( n + 1 \right) \left( n + 2 \right) / 2$ . He did this by finding n + 1 further polynomials that lie in the same space in such a way that the augmented set of polynomials is still linearly independent. More applications of the dimension argument can be found in Graham, Grötschel, and Lovász (1995, chapter 31).

Spectral techniques, that is, an analysis of eigenvectors and eigenvalues [I.3 §4.3], have been used extensively in graph theory. The link comes through the notion of an adjacency matrix of a graph G. This is defined to be the matrix A with entries $_ { a _ { u , \nu } }$ for each pair of (not-necessarily-distinct) vertices u and v, where $a _ { u , \nu } ~ = ~ 1$ if u and v are joined by an edge, and $a _ { u , \nu } = 0$ otherwise. This matrix is symmetric, and therefore, by standard results in linear algebra, it has real eigenvalues and an orthonormal basis [III.37] of eigenvectors. It turns out that there is a tight relationship between the eigenvalues of the adjacency matrix A and several structural properties of the graph G, and these properties can often be useful in the study of various extremal problems. Of particular interest is the second largest eigenvalue of a regular graph. Suppose that every vertex of a graph G has degree d. Then the vector for which every entry is 1 is easily seen to be an eigenvector with eigenvalue d, and this is the largest eigenvalue. If all other eigenvalues have modulus much smaller than d, then it turns out that G behaves in many ways like a random d-regular graph. In particular, the number of edges inside any set of k of the vertices is roughly the same (provided k is not too small) as one would expect with a random graph. It follows easily that any set of vertices that is not too big has many neighbors among the vertices outside that set. Graphs with the latter property are called expanders [III.24] and have numerous applications in theoretical computer science. Constructing such graphs explicitly is not an easy matter and was at one time a major open problem. Now, however, several constructions are known, based on algebraic tools. See chapter 9 of Alon and Spencer (2000), and its references, for more details.

The application of topological methods in the study of combinatorial objects such as partially ordered sets, graphs, and set systems has already become part of the mathematical machinery commonly used in combinatorics. An early example is Lovász’s proof of Kneser’s conjecture, mentioned in section 2.3. Another example is a result of which the following is a representative special case. Suppose you have a piece of string with 10 red beads, 15 blue beads, and 20 yellow beads on it. Then, no matter what order the beads come in, you can cut the string in at most 12 places and place the resulting segments of beaded string into five piles, each of which contains two red beads, three blue beads, and four yellow beads. The number 12 is obtained by multiplying 4, the number of piles minus 1, by 3, the number of colors. The general case of this result was proved by Alon using a generalization of Borsuk’s theorem. Many additional examples of topological proofs appear in Graham, Grötschel, and Lovász (1995, chapter 34).

# 3 Probabilistic Combinatorics

A wonderful development took place in twentieth-century mathematics when it was realized that it is sometimes possible to use probabilistic reasoning to prove mathematical statements that do not have an obvious probabilistic nature. For example, in the first half of the century, Paley, Zygmund, Erd˝os, Turán, Shannon, and others used probabilistic reasoning to obtain striking results in analysis, number theory, combinatorics, and information theory. It soon became clear that the so-called probabilistic method is a very powerful tool for proving results in discrete mathematics. The early results combined combinatorial arguments with fairly elementary probabilistic techniques, but in recent years the method has been greatly developed, and now it often requires one to apply much more sophisticated techniques. A recent text dealing with the subject is Alon and Spencer (2000).

The applications of probabilistic techniques in discrete mathematics were initiated by Paul Erd˝os, who contributed to the development of the method more than anyone else. One can classify them into three groups.

The first deals with the study of certain classes of random combinatorial objects, like random graphs or random matrices. The results here are essentially results in probability theory, although most of them are motivated by problems in combinatorics. A typical problem is the following: if we pick a graph “at random,” what is the probability that it contains a Hamilton cycle?

The second group consists of applications of the following idea. Suppose you want to prove that a combinatorial structure exists with certain properties. Then one possible method is to choose a structure randomly (from a probability distribution that you are free to specify) and estimate the probability that it has the properties you want. If you can show that this probability is greater than 0, then such a structure exists. Surprisingly often it is much easier to prove this than it is to give an example of a structure that works. For instance, is there a graph with large girth (meaning it has no short cycles) and large chromatic number? Even if “large” means “at least 7,” it is very hard to come up with an example of such a graph. But their existence is a fairly easy consequence of the probabilistic method.

The third group of applications is perhaps the most striking of all. There are many examples of statements that appear to be completely deterministic (even when one is used to the idea of using probability to give existence proofs) but that nevertheless yield to probabilistic reasoning. In the remainder of this section we shall briefly describe some typical examples of each of these three kinds of application.

# 3.1 Random Structures

The systematic study of random graphs was initiated by Erd˝os and Rényi in 1960. The most common way of defining a random graph is to fix a probability p and then to join each pair of vertices with an edge with probability p, with all the choices made independently. The resulting graph is denoted $G ( n , p )$ . (Formally speaking, $G ( n , p )$ is not a graph but a probability distribution, but one often talks about it as though it is a graph that has been produced in a random way.) Given any property, such as “contains no triangles,” we can study the probability that $G ( n , p )$ has that property.

A striking discovery of Erd˝os and Rényi was that many properties of graphs “emerge very suddenly.” Some examples are “contains a Hamilton cycle,” “is not planar,” and “is connected.” These properties are all monotone, which means that if a graph G has the property and you add an edge to $G ,$ then the resulting graph still has the property. Let us take one of these properties and define $f ( p )$ to be the probability that the random graph $G ( n , p )$ has it. Because the property is monotone, $f ( p )$ increases as p increases. What Erd˝os and Rényi discovered was that almost all of this increase happens in a very short time. That is, $f ( p )$ is almost 0 for small p and then suddenly changes very rapidly and becomes almost 1.

Perhaps the most famous and illustrative example of this swift change is the sudden appearance of the socalled giant component. Let us look at $G ( n , p )$ when p has the form $c / n .$ . If $c \ < \ 1$ , then with high probability all the connected components of $G ( n , p )$ have size at most logarithmic in n. However, if $c > 1$ , then $G ( n , p )$ ) almost certainly has one component of size linear in n (the giant component), while all the rest have logarithmic size. This is related to the phenomenon of phase transitions in mathematical physics, which are discussed in probabilistic models of critical phenomena [IV.25]. A result of Friedgut shows that the phase transition for a graph property that is “global,” in a sense that can be made precise, is sharper than the one for a “local” property.

Another interesting early discovery in the study of random graphs was that many of the basic parameters of graphs are highly “concentrated.” A striking example that illustrates what this means is the fact that, for any fixed value of p and for most values of $n ,$ almost all graphs $G ( n , p )$ have the same clique number. That is, there exists some r (depending on p and n) such that with high probability, when n is large, the clique number of $G ( n , p )$ is equal to r . Such a result cannot hold for all n, for continuity reasons, but in the exceptional cases there is still some r such that the clique number is almost certainly equal either to r or to r 1. In both cases, r is roughly 2 log n/ log $( 1 / p )$ . The proof of this result is based on the so-called second moment method: one estimates the expectation and the variance of the number of cliques of a given size contained in $G ( n , p )$ , and applies well-known inequalities of Markov and chebyshev [VI.45].

The chromatic number of the random graph $G ( n , p )$ is also highly concentrated. Its typical behavior for values of p that are bounded away from 0 was determined by Bollobás. A more general result, in which p is allowed to tend to 0 as $n  \infty ,$ , was proved by Shamir, Spencer, Łuczak, Alon, and Krivelevich. In particular, it can be shown that for every $\begin{array} { r } { \alpha < \frac { 1 } { 2 } } \end{array}$ and every integervalued function $r ( n ) ~ < ~ n ^ { \alpha }$ , there exists a function $p ( n )$ such that the chromatic number of $G ( n , p ( n ) )$ ) is precisely $r ( n )$ almost surely. However, determining the precise degree of concentration of the chromatic number of $G ( n , p )$ , even in the most basic and important case $\begin{array} { r } { p \ = \ \frac { 1 } { 2 } } \end{array}$ (in which all labeled graphs on n vertices occur with equal probability), remains an intriguing open problem.

Many additional results on random graphs can be found in Janson, Łuczak, and Ruci´nski (2000).

# 3.2 Probabilistic Constructions

One of the first applications of the probabilistic method in combinatorics was a lower bound given by Erd˝os for the Ramsey number $R ( k , k )$ , which was defined in section 2.2. He proved that if

$$
\binom {n} {k} 2 ^ {1 - \binom {k} {2}} <   1,
$$

then $R ( k , k ) > n .$ . That is, there is a red/blue coloring of the edges of the complete graph on n vertices such that no clique of size k is completely red or completely blue. Notice that the number $n = \lfloor 2 ^ { k / 2 } \rfloor$ satisfies the above inequality for all $k \geqslant 3 ,$ so Erd˝os’s result gives an exponential lower bound for $R ( k , k )$ . The proof is simple: if you color the edges randomly and independently, then the probability that any fixed set of k vertices has all its edges of the same color is twice $2 ^ { - { \binom { k } { 2 } } }$ . Thus, the expected number of cliques with this property is

$$
\binom {n} {k} 2 ^ {1 - \binom {k} {2}}.
$$

If this is less than 1, then there must be at least one coloring for which there are no cliques with this property, and the result is proved.

Note that this proof is completely nonconstructive, in the sense that it merely proves the existence of such a coloring, but gives no efficient way of actually constructing one.

A similar computation yields a solution for the tournament problem mentioned in section 1.1. If the results of the tournament are random, then the probability, for any particular k teams, that no other team beats them all is $( 1 - ( 1 / 2 ^ { k } ) ) ^ { n - k }$ . From this it follows that if

$$
\binom {n} {k} \left(1 - \frac {1}{2 ^ {k}}\right) ^ {n - k} <   1,
$$

then there is a nonzero probability that for every choice of k teams, there is another team that beats them all. In particular, it is possible for this to happen. If n is larger than about $k ^ { 2 } 2 ^ { k }$ log 2, then the above inequality holds.

Probabilistic constructions have been very powerful in supplying lower bounds for Ramsey numbers. Besides the bound for $R ( k , k )$ mentioned above, there is a subtle probabilistic proof, due to Kim, that $R ( 3 , k ) \geqslant$ $c k ^ { 2 } /$ log k, for some $c > 0 .$ This is known to be tight up to a constant factor, as proved by Ajtai, Komlós, and Szemerédi, who also used probabilistic methods.

# 3.3 Proving Deterministic Theorems

Suppose that you color the integers with k colors. Let us call a set S multicolored if all k colors appear in S. Straus conjectured that for every k there is an m with the following property: given any set S with m elements, there is a coloring of the integers with k colors such that all translates of S are multicolored. This conjecture was proved by Erd˝os and Lovász. The proof is probabilistic, and applies a tool called the Lovász local lemma, which, unlike many probabilistic techniques, allows one to show that certain events hold with nonzero probability even when this probability is extremely small. The assertion of this lemma, which has numerous additional applications, is, roughly, that for any finite collection of “nearly independent” low-probability events, there is a positive probability that none of the events holds. Note that the statement of Straus’s conjecture has nothing to do with probability, and yet its proof relies on probabilistic arguments.

A graph G is k-colorable, as we have said, if you can properly color its vertices with k colors. Suppose now that instead of trying to use k colors in total, you have a separate list of k colors for each vertex, and this time you want to find a proper coloring of G where each vertex gets a color from its own list. If you can always do so, no matter what the lists are, then G is called k-choosable, and the smallest k for which G is k-choosable is called the choice number ch(G). If all the lists are the same, then one obtains a k-coloring, so ch(G) must be at least as big as $\chi ( G )$ . One might expect ch(G) to be equal to $\chi ( G )$ , since it seems as though using different lists of k colors for different vertices would make it easier to find a proper coloring than using the same k colors for all vertices. However, this turns out to be far from true. It can be proved that for any constant c there is a constant C such that any graph with average degree at least C has choice number at least c. Such a graph might easily be bipartite (and therefore have chromatic number 2), so it follows that ch(G) can be much bigger than $\chi ( G )$ . Somewhat surprisingly, the proof of this result is probabilistic.

An interesting application of this fact concerns a graph that arises in Ramsey theory. Its vertices are all the points in the plane, with two vertices joined by an edge if and only if the distance between them is 1. The choice number of this graph is infinite, by the above result, but the chromatic number is known to be between 4 and 7.

A typical problem in Ramsey theory asks for a substructure of some kind that is entirely colored with one color. Its cousin, discrepancy theory, merely asks that the numbers of times the colors are used are not too close to each other. Probabilistic arguments have proved extremely useful in numerous problems of this general kind. For example, Erd˝os and Spencer proved that in any red/blue coloring of the edges of the complete graph $K _ { n }$ there is a subset $V _ { 0 }$ of vertices such that the difference between the number of red edges inside $V _ { 0 }$ and the number of blue edges inside $V _ { 0 }$ is at least $c n ^ { 3 / 2 }$ , for some absolute constant $c > 0 .$ . This problem is a convincing manifestation of the power of probabilistic methods, since they can be used in the other direction as well, to prove that the result is tight up to a constant factor. Additional examples of such results can be found in Alon and Spencer (2000).

# 4 Algorithmic Aspects and Future Challenges

As we have seen, it is one matter to prove that a certain combinatorial structure exists, and quite another to construct an example. A related question is whether an example can be generated by means of an efficient algorithm [IV.20 §2.3], in which case we call it explicit. This question has become increasingly important because of the rapid development of theoretical computer science, which has close connections with discrete mathematics. It is particularly interesting when the structures in question have been proved to exist by means of probabilistic arguments. Efficient algorithms for producing them are not just interesting on their own, but also have important applications in other areas. For example, explicit constructions of error-correcting codes that are as good as random ones are of major interest in coding and information theory [VII.6], and explicit constructions of certain Ramsey-type colorings may have applications in derandomization [IV.20 §7.1.1] (the process of converting randomized algorithms into deterministic ones).

It turns out, however, that the problem of finding a good explicit construction is often very difficult. Even the simple proof of Erd˝os, described in section 3.2, that there are red/blue colorings of graphs with $\lfloor 2 ^ { k / 2 } \rfloor$ vertices containing no monochromatic clique of size k leads to an open problem that seems very difficult. Can we construct, explicitly, such a graph with $n \geqslant ( 1 + \epsilon ) ^ { k }$ vertices in time that is polynomial in n? Here we allow  to be any constant, as long as it is positive. This problem is still wide open, despite considerable efforts from many mathematicians.

The application of other advanced tools, such as algebraic and analytic techniques, spectral methods, and topological proofs, also tends to lead in many cases to nonconstructive proofs. The conversion of these to algorithmic arguments may well be one of the main future challenges of the area.

Another interesting recent development is the increased appearance of computer-aided proofs in combinatorics, starting with the proof of the four-color theorem [V.12]. To incorporate such proofs into the area, without threatening its special beauty and appeal, is a further challenge.

These challenges, the fundamental nature of the area, its tight connection to other disciplines, and its many fascinating open problems ensure that combinatorics will continue to play an essential role in the general development of mathematics and science in the future.

# Further Reading

Alon, N., and J. H. Spencer. 2000. The Probabilistic Method, 2nd edn. New York: John Wiley.   
Bollobás, B. 1978. Extremal Graph Theory. New York: Academic Press.   
Graham, R. L., M. Grötschel, and L. Lovász, eds. 1995. Handbook of Combinatorics. Amsterdam: North-Holland.   
Graham, R. L., B. L. Rothschild, and J. H. Spencer. 1990. Ramsey Theory, 2nd edn. New York: John Wiley.   
Janson, S., T. Łuczak, and A. Ruci´nski. 2000. Random Graphs. New York: John Wiley.   
Jukna, S. 2001. Extremal Combinatorics. New York: Springer.

Matoušek, J. 2002. Lectures on Discrete Geometry. New York: Springer.   
Nathanson, M. 1996. Additive Number Theory: Inverse Theorems and the Geometry of Sumsets. New York: Springer.   
Pach, J., and P. Agarwal. 1995. Combinatorial Geometry. New York: John Wiley.   
Tao, T., and V. H. Vu. 2006. Additive Combinatorics. Cambridge: Cambridge University Press.

# IV.20 Computational Complexity Oded Goldreich and Avi Wigderson

# 1 Algorithms and Computation

This article is concerned with what can be computed efficiently, and what cannot. We will introduce several important concepts and research areas, such as formal models of computation, measures of efficiency, the P versus NP question, NP-completeness, circuit complexity, proof complexity, randomized computation, pseudorandomness, probabilistic proof systems, cryptography, and more. Underlying them all are the related notions of algorithms and computation, and we begin by discussing these.

# 1.1 What Is an Algorithm?

Suppose that you are presented with a large positive integer N and asked to determine whether it is prime. What should you do? One possibility would be to apply the method of trial division. That is, first see whether N is even, then whether it is a multiple of 3, then whether it is a multiple of 4, and so on through all the numbers up to √N. If N is composite, then it has a factor between 2 and $\sqrt { N } ,$ , so it is prime if and only if the answer to all these questions is no.

The trouble with this method is that it is highly inefficient. Suppose, for instance, that N has 101 digits. Then $\sqrt { N }$ is at least $1 0 ^ { 5 0 }$ , so in order to carry the method out one would have to answer $1 0 ^ { 5 0 }$ questions of the form, “Is K a factor of $N ? ^ { \prime }$ This would take far longer than a human lifetime, even if all the world’s computers devoted themselves to the task. What, then, is an “efficient procedure”? This question divides into two parts: what is a procedure, and what counts as efficient? We shall look at these two questions in turn.

Two very obvious conditions that a method should satisfy if it is to count as a procedure for solving this problem are finiteness—that the procedure should have a finite description (so, for example, one cannot simply look up the answer in an infinite list of integers and their factorizations)—and correctness—that, for every N, it correctly tells you whether N is prime.

There is also a third, more subtle, condition, which goes to the heart of what is meant by the word “algorithm.” It is that it should consist of simple steps. This is needed in order to rule out ridiculous “procedures” such as, “See whether N has any nontrivial factors; declare N to be prime if and only if it does not.” The problem with this is that we cannot see, just like that, whether N has nontrivial factors. By contrast, all that the method of trial division asks of us is that we should do basic arithmetic, such as increasing integers by 1, comparing them, and doing long division. Moreover, the procedures of basic arithmetic can be broken down into yet simpler steps: for instance, it is possible to do long division by a succession of elementary operations applied to single digits at a time.

In order to understand this simplicity condition better, and to prepare ourselves for a formal definition of the notion of algorithms, let us look at long division in slightly more detail. Suppose that you have a piece of paper in front of you and you want to divide 5 959 578 by 857. You will write the two numbers down, and then, as the calculation proceeds, you will write other numbers as well. For instance, you may wish to start by writing out all the multiples of 857 up to 9 × 857. At some point early on you will probably find yourself comparing 5999 = 7 × 857 with 5959: this you do by scanning the numbers from left to right and comparing individual digits. In this case, a difference is first detected in the third digit. You then write 5142 (which is 6  857) underneath the 5959, subtract (again by scanning numbers from left to right and performing single-digit operations), write down the difference 817, “bring down” the next digit, 5, of 5 959 578, and repeat the process with the number 8175.

At each stage in this calculation you are modifying the piece of paper in front of you. As you do so you need to keep track of which stage of the procedure you are at (whether you are writing out the initial table of multiples of 857, or seeing which one is the largest that does not exceed another number, or subtracting one number from another, or bringing down a digit, etc.), and which symbols on the page you are currently dealing with. What is remarkable is that this information has a fixed size, in the sense that it does not increase as the size of the input (that is, the two numbers to be divided) increases.

Therefore, the procedure can be regarded as making local changes to some “environment,” using repeated applications of a fixed rule that does not depend on the input. (This rule will typically have some internal structure, such as a list of simpler rules together with specifications of the circumstances under which they should be applied.) In general, this is what we mean by a computation: it modifies an environment by means of repeated applications of a fixed rule. The rule is usually referred to as an algorithm. Notice that this description applies to many scientific theories of dynamic evolution in nature (of weather, chemical reactions, or biological processes, for example). Thus, these can be regarded as computational processes, of sorts. Some of these dynamical systems also demonstrate well the fact that simple, local rules can result in a very complex modification of the environment if they are iterated many times. (See dynamics [IV.14] for further discussion of this phenomenon.)

Thoughts such as these lie behind the idea of a Turing machine, turing’s [VI.94] famous formalization of the notion of an algorithm. It is interesting that he came up with his formalization before computers existed. Indeed, this abstraction and central features of it, most notably the existence of a “universal” machine, greatly influenced the actual construction of computers.

It is very important to know that the idea of an algorithm can be formalized, so that one can talk precisely about whether there are algorithms that will perform particular tasks, how many steps they need for a given size of input, and so on. However, there are many ways of doing this, which all turn out to be equivalent, and for the purposes of understanding this article it is not necessary to go into the details of any particular method. (You can, if you like, think of an algorithm as any procedure that can be programmed on a real computer—slightly idealized so that it has unlimited storage space—and a step of an algorithm as any change of one of the bits of that computer from a 0 to a 1 or vice versa.) Nevertheless, just to show roughly how it is done, here is a brief description of the basic features of the Turing machine model.

To begin with, one makes the observation that all computational problems can be encoded as operations on sequences of 0s and 1s. (This observation is not just theoretically useful but also very important for the actual building of computers.) For example, all numbers that occur in the course of a computation can be converted into their binary representations; one can also use 1 to stand for “true” and 0 to stand for “false” and thereby perform the basic logical operations; and so on. For this reason we can define a very simple “environment” for a Turing machine: it is a “tape,” infinitely long in both directions, that consists of a row of “cells,” each of which contains either a 0 or a 1. Before the computation starts, a certain prespecified portion of this tape is filled with the input, which is a sequence of 0s and 1s. The algorithm is a little control mechanism. At any one time, this mechanism can be in one of a finite set of states, and it is located at one of the cells of the tape. According to the state it is in and the value, 0 or 1, that it sees at the cell it has reached, it makes three decisions: whether to change the value in the cell, whether to move left or right by one cell, and which state it should next be in.

One of the states of this control mechanism is “halt.” If this state is reached, then the mechanism stops doing anything and is said to have halted. At that point, a certain prespecified portion of the tape will be regarded as the output of the machine. An algorithm can be thought of as any Turing machine that halts for every possible input. And the number of steps of the algorithm is the number of steps taken by that Turing machine. Remarkably, this very simple computational model is enough to capture the full power of computation: in theory one could build a Turing machine, out of clockwork, say, that would be able to do whatever a modern supercomputer can do. (However, it would take too long over each step to be practical for anything but the very simplest of computations.)

# 1.2 What Does an Algorithm Compute?

A Turing machine converts a sequence of 0s and 1s into another sequence of 0s and 1s. If we wish to use mathematical language to discuss this, then we need to give a name to the set of 0, 1 -sequences. To be precise, we consider the set of all finite sequences of 0s and 1s, and we call this set I. It is also useful to write $\mathtt { I } _ { n }$ for the set of all 0, 1 -sequences of length n. If x is a sequence in I, then we write x for its length: for instance, if x is the string 0100101, then $| x | = 7$ . To say that a Turing machine converts a sequence of 0s and 1s into another such sequence (if it halts) is to say that it naturally defines a function from I to I. If M is the Turing machine and fM is the corresponding function, then we say that M computes $f _ { M }$ .

Thus, every function $f : \mathbb { I }  \mathbb { I }$ gives rise to a computational task, namely that of computing f . We say that $f$ is computable if this is possible: that is, if there exists a Turing machine M such that the corresponding function $f _ { M }$ is equal to f . A central early result (due to Turing and independently to church [VI.89]) is that some natural functions are not computable. (For more details, see the insolubility of the halting problem [V.20].) However, complexity theory deals only with computable functions, and studies which of these can be computed efficiently.

Using the notation we have just introduced, we can formally describe various different kinds of computational tasks, of which two major examples are search problems and decision problems. The aim of a search problem is, informally speaking, to find a mathematical object with certain properties: for instance, one might wish to find a solution to a system of equations, and this solution might not be unique. We can model this by means of a binary relation [I.2 §2.3] R on the set I: for a pair $( x , y )$ of strings in I, we say that y is a valid solution of problem instance x if xRy. (This notation means that x is related to y in the way specified by $R ;$ another common notation for the same thing is $( x , y ) \in R . )$ For example, we might let x and y be binary expansions of positive integers N and K, respectively, and say that $x R y$ if and only if N is a composite number and K is a nontrivial factor of N. Informally, this search problem would be, “Find a nontrivial factor of N.” If M is an algorithm that computes a certain function $f _ { M } : \mathbb { I } \to \mathbb { I }$ , then we say that M solves the search problem R if $f _ { M } ( x )$ is a valid solution of x for every problem instance x that has a solution. For example, it solves the search problem just defined if, for every composite number N with binary expansion x, $f _ { M } ( x )$ is the binary expansion of a nontrivial factor K of N.

Notice that in the above example we were interested in positive integers, but formally speaking an algorithm is a function of binary strings. This was not a problem, because there is a convenient and natural way to encode integers as binary strings—via their usual binary expansions. For the rest of this article, we shall feel free to blur the distinction between the mathematical objects we wish to investigate and the strings we use to represent them in a computation. For instance, it is simpler to think of the algorithm M in the previous paragraph as computing a function $f _ { M } : \mathbb { N } \to \mathbb { N } ,$ , and solving the search problem if, for every composite number $N ,$ , $f _ { M } ( N )$ is a nontrivial factor of N. We stress that the representation of objects by strings is a rather succinct one: it takes only log N bits to represent the number $N ,$ so the number N is exponentially larger than the length of its representation.