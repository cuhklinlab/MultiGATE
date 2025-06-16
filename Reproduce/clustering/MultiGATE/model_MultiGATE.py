import tensorflow.compat.v1 as tf
tf.disable_v2_behavior()

class MGATE():

    def __init__(self, hidden_dims1, hidden_dims2, spot_num,  temp = 1.0, nonlinear=True, weight_decay=0.0001):
        self.n_layers = len(hidden_dims1) - 1
     
        self.W1, self.v1, self.prune_v1 = self.define_weights(hidden_dims1, 1)
        self.W2, self.v2, self.prune_v2 = self.define_weights(hidden_dims2, 2)
        self.vgp = self.define_weights_gp(spot_num, 'g')
        self.C1 = {}
        self.prune_C1 = {}
        self.C2 = {}
        self.prune_C2= {}
        self.Cgp = {}
        # self.Cp = {}
        self.nonlinear = nonlinear
        self.weight_decay = weight_decay
        self.hidden_dims1 = hidden_dims1
        self.hidden_dims2 = hidden_dims2
        self.temp = temp

    def __call__(self, A, prune_A, GP, X1, X2):
        # Encoder
        H = tf.concat([tf.transpose(X1), tf.transpose(X2)], axis=0)
        H = self.__encoder_gp(GP, H, self.Cgp, self.vgp, 'gp')
        H = tf.nn.relu(H)
        split_matrices = tf.split(H, num_or_size_splits=[tf.shape(X1)[1], tf.shape(X2)[1]], axis=0) # [tf.shape(X1)[1], tf.shape(X2)[1]]
        H1 = tf.transpose(split_matrices[0])
        H2 = tf.transpose(split_matrices[1]) # X2 #tf.transpose(split_matrices[1]);
        for layer in range(self.n_layers):
            H1 = self.__encoder(A, prune_A, H1, self.W1, self.C1, self.prune_C1, self.v1, self.prune_v1, layer, 1)
            H2 = self.__encoder(A, prune_A, H2, self.W2, self.C2, self.prune_C2, self.v2, self.prune_v2, layer, 2)
            if self.nonlinear:
                if layer != self.n_layers - 1:
                    H1 = tf.nn.elu(H1)
                    H2 = tf.nn.elu(H2)
        # Final node representations
        self.H1 = H1
        self.H2 = H2

        # Decoder
        for layer in range(self.n_layers - 1, -1, -1):
            H1 = self.__decoder(H1, self.W1, self.C1, self.prune_C1, layer)
            H2 = self.__decoder(H2, self.W2, self.C2, self.prune_C2, layer)
            if self.nonlinear:
                # H1 = tf.nn.elu(H1)
                # H2 = tf.nn.elu(H2)
                if layer != 0:
                    H1 = tf.nn.elu(H1)
                    H2 = tf.nn.elu(H2)
        H = tf.concat([tf.transpose(H1), tf.transpose(H2)], axis=0)
        H = self.__decoder_gp(H, self.Cgp)
        H = tf.nn.elu(H)
        split_matrices = tf.split(H, num_or_size_splits=[tf.shape(X1)[1], tf.shape(X2)[1]], axis=0)
        H1 = tf.transpose(split_matrices[0])
        H2 = tf.transpose(split_matrices[1])

        X1_ = H1
        X2_ = H2

        # Define the projection matrices as TensorFlow variables
        W_i = tf.Variable(tf.random_normal(shape=(self.hidden_dims1[2], self.hidden_dims1[2])))
        W_t = tf.Variable(tf.random_normal(shape=(self.hidden_dims2[2], self.hidden_dims2[2])))

        # Define the temperature parameter as a TensorFlow variable
        temperature = tf.Variable(initial_value= self.temp*1.0) #-5.0

        # Define the input tensors for images and text
        # I = tf.placeholder(dtype=tf.float32, shape=(None, ...))
        # T = tf.placeholder(dtype=tf.float32, shape=(None, ...))

        # Extract feature representations of each modality
        # I_f = image_encoder(I)  # [n, d_i]
        # T_f = text_encoder(T)  # [n, d_t]

        # Joint multimodal embedding [n, d_e]
        RNA_e = tf.nn.l2_normalize(tf.matmul(self.H1, W_i), dim=1)
        ATAC_e = tf.nn.l2_normalize(tf.matmul(self.H2, W_t), dim=1)
        self.H1 = RNA_e  # add 1219
        self.H2 = ATAC_e# add 1219

        # Scaled pairwise cosine similarities [n, n]
        logits = tf.matmul(RNA_e, ATAC_e, transpose_b=True) * tf.exp(temperature)
        transposed_logits = tf.transpose(logits)
        # Symmetric loss function
        labels = tf.range(tf.shape(logits)[0])
        loss_rna = tf.nn.sparse_softmax_cross_entropy_with_logits(labels=labels, logits=logits, name="loss_rna")
        loss_atac = tf.nn.sparse_softmax_cross_entropy_with_logits(labels=labels, logits=transposed_logits, name="loss_atac")
        clip_loss = tf.reduce_mean((loss_rna + loss_atac) / 2.0, name="clip_loss")
        #clip_loss = tf.multiply(5.0, clip_loss)  # Multiply clip_loss by 5

        # I_e = l2_normalize(np.dot(self.H1, W_i), axis=1)
        # T_e = l2_normalize(np.dot(self.H2, W_t), axis=1)
        # # scaled pairwise cosine similarities [n, n]
        # logits = np.dot(I_e, T_e.T) * np.exp(t)
        # # symmetric loss function
        # labels = np.arange(n)
        # loss_RNA = cross_entropy_loss(logits, labels, axis=0)
        # loss_ATAC = cross_entropy_loss(logits, labels, axis=1)
        # CLIP_loss = (loss_RNA + loss_ATAC) / 2

        # I_e = l2_normalize(np.dot(I_f, W_i), axis=1)
        # T_e = l2_normalize(np.dot(T_f, W_t), axis=1)
        # # scaled pairwise cosine similarities [n, n]
        # logits = np.dot(I_e, T_e.T) * np.exp(t)
        # # symmetric loss function
        # labels = np.arange(n)
        # loss_i = cross_entropy_loss(logits, labels, axis=0)
        # loss_t = cross_entropy_loss(logits, labels, axis=1)
        # loss = (loss_i + loss_t) / 2

        # The reconstruction loss of node features
        dim_ratio = tf.cast(tf.sqrt(tf.shape(X2)[1]/tf.shape(X1)[1]), dtype=tf.float32) #tf.sqrt()
        features_loss1 = tf.sqrt(tf.reduce_sum(tf.reduce_sum(tf.pow(X1 - X1_, 2))))
        features_loss2 = tf.sqrt(tf.reduce_sum(tf.reduce_sum(tf.pow(X2 - X2_, 2))))
        # features_loss1 = tf.sqrt(tf.reduce_mean(tf.reduce_mean(tf.pow(X1 - X1_, 2))))
        # features_loss2 = tf.sqrt(tf.reduce_mean(tf.reduce_mean(tf.pow(X2 - X2_, 2))))
        

        # for layer in range(self.n_layers):
        #     weight_decay_loss = 0
        #     weight_decay_loss += tf.multiply(tf.nn.l2_loss(self.W[layer]), self.weight_decay, name='weight_loss')

        # Total loss
        # self.loss = features_loss + weight_decay_loss

        # H = X
        # for layer in range(self.n_layers):
        #     H = self.__encoder(A, prune_A, H, layer)
        #     if self.nonlinear:
        #         if layer != self.n_layers - 1:
        #             H = tf.nn.elu(H)
        # # Final node representations
        # self.H = H
        #
        # # Decoder
        # for layer in range(self.n_layers - 1, -1, -1):
        #     H = self.__decoder(H, layer)
        #     if self.nonlinear:
        #         if layer != 0:
        #             H = tf.nn.elu(H)
        # X_ = H
        #
        # # The reconstruction loss of node features
        # features_loss = tf.sqrt(tf.reduce_sum(tf.reduce_sum(tf.pow(X - X_, 2))))
        #
        weight_decay_loss = 0
        for layer in range(self.n_layers):
            
            weight_decay_loss += tf.multiply(tf.nn.l2_loss(self.W1[layer]), self.weight_decay, name='weight_loss')
            weight_decay_loss += tf.multiply(tf.nn.l2_loss(self.W2[layer]), self.weight_decay, name='weight_loss')
        weight_decay_loss += tf.multiply(tf.nn.l2_loss(W_i), self.weight_decay, name='weight_loss')
        weight_decay_loss += tf.multiply(tf.nn.l2_loss(W_t), self.weight_decay, name='weight_loss')

        #
        # # Total loss
        # features_loss1*10
        # features_loss1=features_loss1*10
        # features_loss2=features_loss2*10
        self.loss = features_loss1   + features_loss2 + weight_decay_loss + clip_loss # *0.0
        self.loss_rna = features_loss1
        self.loss_atac = features_loss2
        self.weight_decay_loss = weight_decay_loss
        self.clip_loss = clip_loss
        # print("features_loss_rna:,  %.4f" % features_loss1.numpy()[0])
        # print("features_loss_atac:,  %.4f" % features_loss2.numpy()[0])
        # print("weight_decay_loss:,  %.4f" % weight_decay_loss.numpy()[0])
        # print("clip_loss:,  %.4f" % clip_loss.numpy()[0])

       
        # self.Att_l = self.C1 # {'C1': self.C1, 'C2': self.C2}
        # self.Att_l = {'C1': self.C1, 'prune_C1': self.prune_C1, 'C2': self.C2, 'prune_C2': self.prune_C2}
        self.Att_l1 = self.C1
        self.Att_l2 = self.C2
        self.Att_lgp = self.Cgp
        
        return self.loss, self.loss_rna, self.loss_atac, self.weight_decay_loss, self.clip_loss, self.H1, self.H2, \
               self.Att_l1, self.Att_l2, self.Att_lgp, X1_, X2_
    def __encoder_gp(self, A, H, C, v, k):
        # H = tf.matmul(H, W)
        C[0] = self.graph_attention_layer(A, H, v, 0, k)
        return tf.sparse_tensor_dense_matmul(C[0], H)


    def __encoder(self, A, prune_A, H, W, C, prune_C, v, prune_v, layer, k):
        H = tf.matmul(H, W[layer])
        if layer == self.n_layers - 1:
            return H
        C[layer] = self.graph_attention_layer(A, H, v[layer], layer, k)

        return tf.sparse_tensor_dense_matmul(C[layer], H)


    def __decoder_gp(self,  H,  C):
        # H = tf.matmul(H, W, transpose_b=True)
        return  tf.sparse_tensor_dense_matmul(C[0], H)

    def __decoder(self, H, W, C, prune_C, layer):
        H = tf.matmul(H, W[layer], transpose_b=True)
        if layer == 0:
            return H

        return tf.sparse_tensor_dense_matmul(C[layer - 1], H)
       

    # def __encoder1(self, A, prune_A, H, layer):
    #     H = tf.matmul(H, self.W[layer])
    #     if layer == self.n_layers - 1:
    #         return H
    #     self.C[layer] = self.graph_attention_layer(A, H, self.v[layer], layer)
    #     if self.alpha == 0:
    #         return tf.sparse_tensor_dense_matmul(self.C[layer], H)
    #     else:
    #         self.prune_C[layer] = self.graph_attention_layer(prune_A, H, self.prune_v[layer], layer)
    #         return (1 - self.alpha) * tf.sparse_tensor_dense_matmul(self.C[layer],
    #                                                                 H) + self.alpha * tf.sparse_tensor_dense_matmul(
    #             self.prune_C[layer], H)
    #
    # def __decoder1(self, H, layer):
    #     H = tf.matmul(H, self.W[layer], transpose_b=True)
    #     if layer == 0:
    #         return H
    #     if self.alpha == 0:
    #         return tf.sparse_tensor_dense_matmul(self.C[layer - 1], H)
    #     else:
    #         return (1 - self.alpha) * tf.sparse_tensor_dense_matmul(self.C[layer - 1],
    #                                                                 H) + self.alpha * tf.sparse_tensor_dense_matmul(
    #             self.prune_C[layer - 1], H)
    #
    # def __encoder2(self, A, prune_A, H, layer):
    #     H = tf.matmul(H, self.W[layer])
    #     if layer == self.n_layers - 1:
    #         return H
    #     self.C[layer] = self.graph_attention_layer(A, H, self.v[layer], layer)
    #     if self.alpha == 0:
    #         return tf.sparse_tensor_dense_matmul(self.C[layer], H)
    #     else:
    #         self.prune_C[layer] = self.graph_attention_layer(prune_A, H, self.prune_v[layer], layer)
    #         return (1 - self.alpha) * tf.sparse_tensor_dense_matmul(self.C[layer],
    #                                                                 H) + self.alpha * tf.sparse_tensor_dense_matmul(
    #             self.prune_C[layer], H)
    #
    # def __decoder2(self, H, layer):
    #     H = tf.matmul(H, self.W[layer], transpose_b=True)
    #     if layer == 0:
    #         return H
    #     if self.alpha == 0:
    #         return tf.sparse_tensor_dense_matmul(self.C[layer - 1], H)
    #     else:
    #         return (1 - self.alpha) * tf.sparse_tensor_dense_matmul(self.C[layer - 1],
    #                                                                 H) + self.alpha * tf.sparse_tensor_dense_matmul(
    #             self.prune_C[layer - 1], H)
    def define_weights_gp(self, spot_num, k):
        # W = tf.get_variable("W%s" % 0 + str(k), shape=(spot_num, spot_num)) #  var_name

        v = {}
        v[0] = tf.get_variable("v%s_0" % 0 + str(k), shape=(spot_num, 1)) #   var_name+"_0"
        v[1] = tf.get_variable("v%s_1" % 0 + str(k), shape=(spot_num, 1)) #  var_name+"_1"

        
        return v

    def define_weights(self, hidden_dims, k):
        W = {}

        for i in range(self.n_layers):
            # var_name = "W{}_{}".format(k, i)
            W[i] = tf.get_variable("W%s" % i + str(k), shape=(hidden_dims[i], hidden_dims[i + 1])) #  var_name

        Ws_att = {}
        for i in range(self.n_layers - 1):
            # var_name = "v{}_{}".format(k, i)
            v = {}
            v[0] = tf.get_variable("v%s_0" % i + str(k), shape=(hidden_dims[i + 1], 1)) #   var_name+"_0"
            v[1] = tf.get_variable("v%s_1" % i + str(k), shape=(hidden_dims[i + 1], 1)) #  var_name+"_1"

            Ws_att[i] = v
        
        return W, Ws_att, None
        prune_Ws_att = {}
        for i in range(self.n_layers - 1):
            # var_name = "prune_v{}_{}".format(k, i)
            prune_v = {}
            prune_v[0] = tf.get_variable("prune_v%s_0" % i + str(k), shape=(hidden_dims[i + 1], 1)) # "prune_v%s_0" % i  var_name+"_0"
            prune_v[1] = tf.get_variable("prune_v%s_1" % i + str(k), shape=(hidden_dims[i + 1], 1)) # "prune_v%s_1" % i  var_name+"_1"

            prune_Ws_att[i] = prune_v

        return W, Ws_att, prune_Ws_att

    def graph_attention_layer(self, A, M, v, layer, k):
        B = tf.SparseTensor(indices=A.indices, values=tf.ones(tf.shape(A.indices)[0]), dense_shape=A.dense_shape)        
        with tf.variable_scope("layer_%s" % layer + str(k)):
            f1 = tf.matmul(M, v[0])
            f1 = B * f1
            f2 = tf.matmul(M, v[1])
            f2 = B * tf.transpose(f2, [1, 0])
            logits = tf.sparse_add(f1, f2)

            unnormalized_attentions = tf.SparseTensor(indices=logits.indices,
                                                      values=  tf.log(A.values) * tf.nn.sigmoid(logits.values), #tf.log(A.values) *
                                                      dense_shape=logits.dense_shape)
            #exponential_weights = tf.SparseTensor(indices=A.indices,
            #                                           values=tf.exp(A.values),
            #                                           dense_shape=A.dense_shape)
            attentions = tf.sparse_softmax(unnormalized_attentions)

            attentions = tf.SparseTensor(indices=attentions.indices,
                                         values= attentions.values, #A.values *
                                         dense_shape=attentions.dense_shape)

            return attentions