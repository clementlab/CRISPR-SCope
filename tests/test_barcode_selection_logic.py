# tests/test_selection_logic.py
import pandas as pd

def test_selection_via_isin_behaves():
    # simulate parsed_information produced by stratify_data
    df = pd.DataFrame({
        'Barcode': ['A','B','C','D'],
        'Color': ['HQ_HI','LQ_HI','HQ_LO','LQ_LO']
    }).set_index('Barcode')

    # single selection
    sel = ['HQ_HI']
    picked = df.loc[df['Color'].isin(sel)].index.tolist()
    assert picked == ['A']

    # multi selection
    sel = ['HQ_HI','HQ_LO']
    picked = set(df.loc[df['Color'].isin(sel)].index.tolist())
    assert picked == {'A','C'}

    # empty selection
    sel = []
    picked = df.loc[df['Color'].isin(sel)].index.tolist()
    assert picked == []

