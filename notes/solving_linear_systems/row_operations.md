--
layout: page
title: Row Operations
nav_order: 7
has_children: false
has_toc: false
parent: Solving Linear Systems
grand_parent: Notes
---

## Changing the system in controlled ways

We have already seen [row echelon form]({% link notes/solving_linear_systems/row_echelon_form.md %}). Our goal now is to identify transformations 
of the system, or equivalently the augmented matrix that reduces the coefficient matrix to 
row echelon form. 

### Scaling a row 

One of the simplest operation is multiplying a row by a _non-zero_ scalar. For example,
let's take the system  
$$
	x + y = 1 \\
	x - y = 3 \\
	2x + y = 0 
$$
and scale the second row by $1/5$. The resulting system 
$$
	x + y = 1 \\
    {} \\
	\frac{1}{5}x - \frac{1}{5}y = \frac{3}{5} \\
    {} \\
	2x + y = 0 
$$
has the same set of solutions as the original system. 

Let's try to capture the idea of "the new system obtained from the old system by 
scaling the a row by a constant" in terms of matrix notation. 

By multiplying by specific matrices, we can achieve a good number of interesting 
operations. 

**Definition**: Let $S(i,c)\_m$ denote the $m \times m$ matrix whose entries are $0$ off the 
diagonal and all of whose entries on the diagonal are $1$ _except_ for the i-th 
entry which is the scalar $c$. Here are a couple examples to make the notation 
clear. 

$$
	S(1,\pi)_2 = \begin{pmatrix} \pi & 0 \\ 0 & 1 \end{pmatrix} \ , \
	S(2,5)_3 = 
	\begin{pmatrix}
		1 & 0 & 0 \\
		0 & 5 & 0 \\
		0 & 0 & 1 
	\end{pmatrix}
$$

As part of [Homework 1]({% link homework/01.md %}), 
you will need to verify the following fact: 

**Lemma**: For any $m \times n$ matrix $A$, the product $B := S(l,c)_m A$ is a 
$m \times n$ matrix satisfying 
$$
	b_{ij} = 
				\begin{cases} 
					a_{ij} & i \neq l \\
					c_{lj} & i = l 
				\end{cases}
$$
In other words, $S(l,c)_m$ achieves through matrix multiplication exactly what we 
want: it scales the l-th row by $c$ and leaves the others unchanged. 

### Exchanging rows

Another simple operation that obviously is swapping the 
two rows (or two equations) in the system. 

Like scaling a row, we can also express exchanging rows through multiplication with a 
special matrix. 

**Definition** The $m \times m$ $(i,j)$ permutation matrix $P(i,j)_m$ is the matrix 
obtained from the $m \times m$ identity by exchanging the row i and row j. 

Some examples 
$$
	P(1,2)_2 = 
					\begin{pmatrix}
						0 & 1 \\
						1 & 0 
					\end{pmatrix}
$$
and 
$$
	P(2,4)_4 = 
					\begin{pmatrix}
						1 & 0 & 0 & 0 \\
						0 & 0 & 0 & 1 \\
						0 & 0 & 1 & 0 \\
						0 & 1 & 0 & 0 
					\end{pmatrix}
$$

**Proposition** For any $m \times n$ matrix $B$ the matrix $P(i,j)_m B$ is obtained 
from $B$ by exchanging row i and row j of $B$. 

**Proof** Let's use some facts about matrix multiplication and clever rewriting of 
$P(i,j)\_m$. We can write 
$$
	P(i,j)_m = I_m - D(i,i) - D(j,j)+ D(j,i) + D(i,j) 
$$
where $I\_m$ is the $m \times m$ identity matrix and $D(i,j)$ is matrix described in 
(broken link)

So 
$$
	P(i,j) B = B - D(i,i) B - D(j,j) B + D(i,j) B + D(i,j) B
$$

From the [Homework 1]({% link homework/01.md %}), we know that $D(i,j) B$ 
is the matrix whose which is $0$ 
except for the i-th row which is equal to the j-th row of $B$.  So reading 
$$
	B - D(i,i) B - D(j,j) B + D(j,i) B + D(i,j) B
$$
we have $B$ then we subtract the delete the i-th and j-th rows 
$$
	B - D(i,i) B - D(j,j) B.
$$
After we add back in the i-th row in the j-th place and the j-th row in the i-th spot,
$$
	B - D(i,i) B - D(j,j) B + D(j,i) B + D(i,j) B.
$$
The result is that we have replaced the i-th row in $B$ with the j-th and also the 
j-th row with the i-th as we wanted to check. 

Thus, we can express the row-exchanged system as the matrix product 
$ P(i,j) (A \mid b)$. 

### Adding multiples of rows to other rows 

For the final operation, it is perhaps less obvious that we preserve the solution set 
of the linear system. 

With permutations, we can make sure all the zero rows are moved to the bottom of 
the matrix. 

With scaling, we can change the pivot to have value $1$. 

Now, we want to use this to zero out the other entries either below (for row echelon)
or both above and below (for reduced row echelon form). 

We do this by adding scalar multiples of one row to another. Let's see this in an 
example.

Consider the matrix
$$
	\begin{pmatrix}
		1 & 2 & 3 \\
		4 & 5 & 6 \\
		7 & 8 & 9
	\end{pmatrix}
$$
Let's subtract $4 \mathbf{R}\_1$ from $\mathbf{R}\_2$. We have 
$$
	4 \mathbf{R}_1 = \begin{pmatrix} 4 & 8 & 12 \end{pmatrix} 
$$
so 
$$
	\mathbf{R}_2 - 4 \mathbf{R}_1 = \begin{pmatrix} 0 & -3 & -6 \end{pmatrix}
$$
Our new matrix is
$$
	\begin{pmatrix}
		1 & 2 & 3 \\
		0 & -3 & -6 \\
		7 & 8 & 9
	\end{pmatrix}
$$

As with scaling a row and switching rows, we want to express the operation $\mathbf{R}_i \mapsto 
\mathbf{R}_i - c\mathbf{R}_j$ via multiplication with a particular matrix. 

Consider the matrix 
$$
	I_m - c D(i,j) 
$$
If we multiply this by a matrix $A$, 
$$
	\left( I_m - c D(i,j) \right) A = I_m A - c D(i,j) A = A - c D(i,j) A
$$
then we are subtracting off a copy of the i-th row in the j-th spot which is scaled by $c$. 

This is exactly what we wanted to do! Thus, the multiplication with $I - c D(i,j)$ subtracts 
$c \mathbf{R}_j^A$ from the $\mathbf{R}_i^A$. 

Let's give these matrices their own names since they will be important:
$$
	L(i,j,c)_m := I_m + cD(i,j)_m 
$$
So for example $L(1,2,10)\_3$ is 
$$
	\begin{pmatrix} 
		1 & 10 & 0 \\
		0 & 1 & 0 \\
		0 & 0 & 1
	\end{pmatrix}
$$

