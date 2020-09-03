from sklearn.feature_extraction.text import TfidfVectorizer  

def top_x_keywords(x, text,**kwargs):
    
    """
    Extract keywords by TFIDF metric
    Keywords for each document (row) are the top x as sorted by TFIDF
    
    Can pass any argument to vectorizer based on the task at hand
    Implented using explicit looping through sparse matrix to avoid memory errors by
    forcing dense representation through numpys sort function
    
    """

    vect = TfidfVectorizer(**kwargs).fit(text)
    TFID_key_words = vect.transform(text)
    tfidf = TFID_key_words.tocoo()
    vocab = vect.get_feature_names()
    top_words = [[] for _ in range(len(text))]

    prev_row = 0
    cx = TFID_key_words.tocoo()
    
    for row,col,val in zip(cx.row, cx.col, cx.data):

        if row != prev_row:

            top_words[prev_row].sort(key = lambda k: k[1], reverse= True) 
            to_check = min(x, len(top_words[prev_row]))
            top_words[prev_row] = [ vocab[top_words[prev_row][i][0]] for i in range(to_check)]

        else:

            top_words[row].append((col,val))


        prev_row = row

    top_words[prev_row].sort(key = lambda k: k[1], reverse= True) 
    to_check = min(x, len(top_words[prev_row]))
    top_words[prev_row] = [ vocab[top_words[prev_row][i][0]] for i in range(to_check)]
    
    return top_words
