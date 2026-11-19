import numpy as np
from numpy import ndarray

def calc_cyclic_matrix_af3(residue_index,
                           head_to_tail: bool = False,
                           list_cid_ss: list[list[int]] = [],
                           token_wise=True,
                           signal_bug_fix=True):

  if list_cid_ss is None:
    list_cid_ss = []
    
  if not head_to_tail and (list_cid_ss is None or list_cid_ss == [] or list_cid_ss == [[]]):
    return [[0]]
  residue_index = np.asarray(residue_index, dtype=np.int32)
  if residue_index.ndim != 1 or residue_index.size == 0:
    raise ValueError('residue_index must be a non-empty one-dimensional array')

  # token_offset is consumed only for pairs of tokens belonging to the same
  # residue.  It must remain the ordinary token-index difference, especially
  # for non-standard residues represented by multiple atom tokens.  Applying
  # residue graph distances here collapses all of those offsets to zero.
  if token_wise:
    token_index = np.arange(residue_index.size, dtype=np.int32)
    return token_index[:, None] - token_index[None, :]

  if len(list_cid_ss) > 1:
    raise ValueError('Only one cyclic peptide chain is currently supported')

  # Non-standard polymer residues can be represented by several atom tokens,
  # all carrying the same residue ID.  Counting occurrences of residue ID 1
  # therefore mistakes the second atom token for the start of another chain.
  # Build the topology in residue space and map it back to token space.
  chain_start = 0
  if list_cid_ss:
    chain_number = list_cid_ss[0][0]
    chain_starts = np.concatenate(
        ([0], np.flatnonzero(np.diff(residue_index) < 0) + 1)
    )
    if not 1 <= chain_number <= len(chain_starts):
      raise ValueError(f'Invalid cyclic peptide chain index: {chain_number}')
    chain_start = int(chain_starts[chain_number - 1])
  chain_end_candidates = np.flatnonzero(
      np.diff(residue_index[chain_start:]) < 0
  )
  chain_end = (
      chain_start + int(chain_end_candidates[0]) + 1
      if chain_end_candidates.size
      else residue_index.size
  )
  chain_residue_index = residue_index[chain_start:chain_end]
  residue_ids = np.unique(chain_residue_index)
  residue_to_offset = {int(res_id): i for i, res_id in enumerate(residue_ids)}

  c1 = []
  c2 = []
  if list_cid_ss:
    topology = list_cid_ss[0]
    if len(topology) < 3 or (len(topology) - 1) % 2:
      raise ValueError('Disulfide topology must be [chain, res1, res2, ...]')
    try:
      c1 = [residue_to_offset[ind] for ind in topology[1::2]]
      c2 = [residue_to_offset[ind] for ind in topology[2::2]]
    except KeyError as exc:
      raise ValueError(
          f'Disulfide residue {exc.args[0]} is absent from chain {topology[0]}'
      ) from exc

  residue_offset = calc_cyclic_matrix_monomer(
      len(residue_ids), c1, c2, head_to_tail=head_to_tail
  )
  token_to_residue = np.array(
      [residue_to_offset[int(res_id)] for res_id in chain_residue_index]
  )
  return residue_offset[np.ix_(token_to_residue, token_to_residue)]

def calc_cyclic_matrix_monomer(n_aa, c1, c2, head_to_tail=False):
  """
  :param n_aa: number of amino acid residues
  :param c1: list of Cys residues
  :param c2: list of the corresponding Cys residues to c1, the same length with c1
  :param head_to_tail:
  :return:
  """
  if len(c1) != len(c2):
    return []  
  # init adjacency matrix
  matrix = np.zeros((n_aa, n_aa)) + n_aa
  for i in range(n_aa):
    matrix[i][i] = 0  
  # linear peptide connection
  for i in range(n_aa - 1):
    matrix[i][i + 1] = 1
    matrix[i + 1][i] = 1  
  # nc connection
  if head_to_tail:
    matrix[0][n_aa - 1] = 1
    matrix[n_aa - 1][0] = 1  
  # ss connection
  for i in range(len(c1)):
    matrix[c1[i]][c2[i]] = 1
    matrix[c2[i]][c1[i]] = 1  
  # get the shortest path
  matrix = calc_shortest_path(matrix)  
  for i in range(matrix.shape[0]):
    for j in range(i + 1, matrix.shape[0], 1):
      matrix[i][j] *= -1  
  return matrix


def calc_shortest_path(matrix):
  """
  Floyd algorithm to find the shortest path
  :param matrix:
  :return:
  """
  path = np.zeros_like(matrix)
  for i in range(matrix.shape[0]):
    path[i] = [j for j in range(matrix.shape[0])]

  for m in range(matrix.shape[0]):
    for i in range(matrix.shape[0]):
      for j in range(matrix.shape[0]):
        if matrix[i][m] + matrix[m][j] < matrix[i][j]:
          matrix[i][j] = matrix[i][m] + matrix[m][j]
          path[i][j] = m
  return matrix


def construct_cyclic_offset_matrix_v20250306(n_res, signal_bug_fix=True,):
  i = np.arange(n_res)
  offset = i[:, None] - i[None, :]
  ij = np.stack([i, i + n_res], -1)
  c_offset = np.abs(ij[:, None, :, None] - ij[None, :, None, :]).min((2, 3))
  if signal_bug_fix:
    a = c_offset < np.abs(offset)
    c_offset[a] = -c_offset[a]
  final_offset = c_offset * np.sign(offset)
  return final_offset
