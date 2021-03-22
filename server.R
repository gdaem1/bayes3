library(shiny)
library(plotly)
source('utilities2.R')

shinyServer(function(input, output) {
  
  observeEvent(input$button, {
    req(
      input$success_A,
      input$total_A,
      input$success_B,
      input$total_B,
      input$success_C,
      input$total_C,
      input$rev_A,
      input$rev_B,
      input$rev_C,
      input$retention_A,
      input$retention_B,
      input$retention_C,
      input$sim_sample
    )
    if(
      input$success_A >= 0 &&
      input$total_A > 0 &&
      input$success_B >= 0 &&
      input$total_B > 0 &&
      input$success_C >= 0 &&
      input$total_C > 0 &&
      input$success_A <= input$total_A &&
      input$success_B <= input$total_B &&
      input$success_C <= input$total_C &&
      input$rev_A >= 0 &&
      input$rev_B >= 0 &&
      input$rev_C >= 0 &&
      input$retention_A >= 0 &&
      input$retention_B >= 0 &&
      input$retention_C >= 0 &&
      input$retention_A <= input$total_A &&
      input$retention_B <= input$total_B &&
      input$retention_C <= input$total_C &&
      input$sim_sample >= 2
    ) {
      sample_A <- isolate({input$total_A})
      sample_B <- isolate({input$total_B})
      sample_C <- isolate({input$total_C})
      conv_A <- isolate({input$success_A/input$total_A})
      conv_B <- isolate({input$success_B/input$total_B})
      conv_C <- isolate({input$success_C/input$total_C})
      arppu_A <- isolate({input$rev_A/input$success_A})
      arppu_B <- isolate({input$rev_B/input$success_B})
      arppu_C <- isolate({input$rev_C/input$success_C})
      arpu_A <- isolate({input$rev_A/input$total_A})
      arpu_B <- isolate({input$rev_B/input$total_B})
      arpu_C <- isolate({input$rev_C/input$total_C})
      alpha_A <- isolate({input$success_A + 1})
      alpha_B <- isolate({input$success_B + 1})
      alpha_C <- isolate({input$success_C + 1})
      beta_A <- isolate({input$total_A - input$success_A + 1})
      beta_B <- isolate({input$total_B - input$success_B + 1})
      beta_C <- isolate({input$total_C - input$success_C + 1})
      k_A <- isolate({input$success_A + 1})
      k_B <- isolate({input$success_B + 1})
      k_C <- isolate({input$success_C + 1})
      theta_A <- isolate({1/(1 + input$rev_A)})
      theta_B <- isolate({1/(1 + input$rev_B)})
      theta_C <- isolate({1/(1 + input$rev_C)})
      convret_A <- isolate({input$retention_A/input$total_A})
      convret_B <- isolate({input$retention_B/input$total_B})
      convret_C <- isolate({input$retention_C/input$total_C})
      alpharet_A <- isolate({input$retention_A + 1})
      alpharet_B <- isolate({input$retention_B + 1})
      alpharet_C <- isolate({input$retention_C + 1})
      betaret_A <- isolate({input$total_A - input$retention_A + 1})
      betaret_B <- isolate({input$total_B - input$retention_B + 1})
      betaret_C <- isolate({input$total_C - input$retention_C + 1})
      retention_A <- isolate({input$retention_A})
      retention_B <- isolate({input$retention_B})
      retention_C <- isolate({input$retention_C})
      res <- isolate({
        bayes_arpu(
          alphaA = alpha_A, betaA = beta_A,
          kA = k_A, thetaA = theta_A,
          alphaB = alpha_B, betaB = beta_B,
          kB = k_B, thetaB = theta_B,
          alphaC = alpha_C, betaC = beta_C,
          kC = k_C, thetaC = theta_C,
          MSamples = input$sim_sample
        )
      })
      post_sample_A <- res$sampleLambdaA/res$sampleOmegaA
      post_sample_B <- res$sampleLambdaB/res$sampleOmegaB
      post_sample_C <- res$sampleLambdaC/res$sampleOmegaC
      diff_post_sampleAB <- post_sample_B - post_sample_A
      diff_post_sampleAC <- post_sample_C - post_sample_A
      diff_post_sampleBC <- post_sample_C - post_sample_B
      hdi_A <- hdi_of_sample(post_sample_A)
      hdi_B <- hdi_of_sample(post_sample_B)
      hdi_C <- hdi_of_sample(post_sample_C)
      hdi_diffAB <- hdi_of_sample(diff_post_sampleAB)
      hdi_diffAC <- hdi_of_sample(diff_post_sampleAC)
      hdi_diffBC <- hdi_of_sample(diff_post_sampleBC)
      x_lim <- {
        a <- min(hdi_A, hdi_B)
        b <- max(hdi_A, hdi_B)
        c(1.2*a - 0.2*b, 1.2*b - 0.2*a)
      }
      x_lim_diff <- {
        a <- hdi_diffAB[1]
        b <- hdi_diffAB[2]
        c(1.2*a - 0.2*b, 1.2*b - 0.2*a)
      }
      x_limAC <- {
        a <- min(hdi_A, hdi_C)
        b <- max(hdi_A, hdi_C)
        c(1.2*a - 0.2*b, 1.2*b - 0.2*a)
      }
      x_lim_diffAC <- {
        a <- hdi_diffAC[1]
        b <- hdi_diffAC[2]
        c(1.2*a - 0.2*b, 1.2*b - 0.2*a)
      }
      x_limBC <- {
        a <- min(hdi_B, hdi_C)
        b <- max(hdi_B, hdi_C)
        c(1.2*a - 0.2*b, 1.2*b - 0.2*a)
      }
      x_lim_diffBC <- {
        a <- hdi_diffBC[1]
        b <- hdi_diffBC[2]
        c(1.2*a - 0.2*b, 1.2*b - 0.2*a)
      }
      
      printPlot <- isolate({TRUE})
    } else {
      sample_A <- isolate({0})
      sample_B <- isolate({0})
      sample_C <- isolate({0})
      conv_A <- isolate({NaN})
      conv_B <- isolate({NaN})
      conv_C <- isolate({NaN})
      arppu_A <- isolate({NaN})
      arppu_B <- isolate({NaN})
      arppu_C <- isolate({NaN})
      arpu_A <- isolate({NaN})
      arpu_B <- isolate({NaN})
      arpu_C <- isolate({NaN})
      alpha_A <- isolate({1})
      alpha_B <- isolate({1})
      alpha_C <- isolate({1})
      beta_A <- isolate({1})
      beta_B <- isolate({1})
      beta_C <- isolate({1})
      k_A <- isolate({1})
      k_B <- isolate({1})
      k_C <- isolate({1})
      theta_A <- isolate({1})
      theta_B <- isolate({1})
      theta_C <- isolate({1})
      hdi_A <- isolate({c(NaN, NaN)})
      hdi_B <- isolate({c(NaN, NaN)})
      hdi_C <- isolate({c(NaN, NaN)})
      hdi_diffAB <- isolate({c(NaN, NaN)})
      hdi_diffAC <- isolate({c(NaN, NaN)})
      hdi_diffBC <- isolate({c(NaN, NaN)})
      x_lim <- isolate({c(NaN, NaN)})
      x_lim_diff <- isolate({c(NaN, NaN)})
      x_limAC <- isolate({c(NaN, NaN)})
      x_lim_diffAC <- isolate({c(NaN, NaN)})
      x_limAD <- isolate({c(NaN, NaN)})
      x_lim_diffBC <- isolate({c(NaN, NaN)})
      x_limBD <- isolate({c(NaN, NaN)})
      
      res <- isolate({
        list(
          convProbBbeatsA = NaN, convExpLossA_AB = NaN, convExpLossB_AB = NaN,
          revProbBbeatsA = NaN, revExpLossA_AB = NaN, revExpLossB_AB = NaN,
          arpuProbBbeatsA = NaN, arpuExpLossA_AB = NaN, arpuExpLossB_AB = NaN,
          convProbAbeatsB = NaN, convExpLossA_AB2 = NaN, convExpLossB_AB2 = NaN,
          revProbAbeatsB = NaN, revExpLossA_AB2 = NaN, revExpLossB_AB2 = NaN,
          arpuProbAbeatsB = NaN, arpuExpLossA_AB2 = NaN, arpuExpLossB_AB2 = NaN,
          
          convProbCbeatsA = NaN, convExpLossA_AC = NaN, convExpLossC_AC = NaN,
          revProbCbeatsA = NaN, revExpLossA_AC = NaN, revExpLossC_AC = NaN,
          arpuProbCbeatsA = NaN, arpuExpLossA_AC = NaN, arpuExpLossC_AC = NaN,
          convProbAbeatsC = NaN, convExpLossA_AC2 = NaN, convExpLossC_AC2 = NaN,
          revProbAbeatsC = NaN, revExpLossA_AC2 = NaN, revExpLossC_AC2 = NaN,
          arpuProbAbeatsC = NaN, arpuExpLossA_AC2 = NaN, arpuExpLossC_AC2 = NaN,
          
          convProbCbeatsB = NaN, convExpLossB_BC = NaN, convExpLossB_BC = NaN,
          revProbCbeatsB = NaN, revExpLossB_BC = NaN, revExpLossC_BC = NaN,
          arpuProbCbeatsB = NaN, arpuExpLossB_BC = NaN, arpuExpLossC_BC = NaN,
          convProbBbeatsC = NaN, convExpLossB_BC2 = NaN, convExpLossB_BC2 = NaN,
          revProbBbeatsC = NaN, revExpLossB_BC2 = NaN, revExpLossC_BC2 = NaN,
          arpuProbBbeatsC = NaN, arpuExpLossB_BC2 = NaN, arpuExpLossC_BC2 = NaN,
        )
      })
      printPlot <- isolate({FALSE})
    }
    
    output$table1 <- renderTable({
      tab <- data.frame(
        metric = c(
          'Sample Size', 'Retained Players','<strong>Conversion<strong>','<strong>Retention<strong>', 'ARPPU',
          'ARPU', '95% HDI'
        ),
        A = c(
          sprintf('\n%.d', sample_A),
          sprintf('\n%.d', retention_A),
          sprintf('\n%.2g%%', conv_A*100),
          sprintf('\n%.2g%%', convret_A*100),
          sprintf('\n%.2g â¬', arppu_A),
          sprintf('\n%.2g â¬', arpu_A),
          sprintf('[%.2g â¬, \n%.2g â¬]', hdi_A[1], hdi_A[2])
        ),
        B = c(
          sprintf('\n%.d', sample_B),
          sprintf('\n%.d', retention_B),
          sprintf('\n%.2g%%', conv_B*100),
          sprintf('\n%.2g%%', convret_B*100),
          sprintf('\n%.2g â¬', arppu_B),
          sprintf('\n%.2g â¬', arpu_B),
          sprintf('[%.2g â¬, \n%.2g â¬]', hdi_B[1], hdi_B[2])
        ),
        C = c(
          sprintf('\n%.d', sample_C),
          sprintf('\n%.d', retention_C),
          sprintf('\n%.2g%%', conv_C*100),
          sprintf('\n%.2g%%', convret_C*100),
          sprintf('\n%.2g â¬', arppu_C),
          sprintf('\n%.2g â¬', arpu_C),
          sprintf('[%.2g â¬, \n%.2g â¬]', hdi_C[1], hdi_C[2])
        )
      )
      colnames(tab) <- c(' ', 'A', 'B', 'C')
      tab
    }, spacing = 'xs', sanitize.text.function = function(x){x})
    
    output$table2 <- renderTable({
      tab <- data.frame(
        column1 = c(
          'Probability that B is better than A',
          'Probability that C is better than A',
          'Probability that C is better than B'
        ),
        conversion = c(
          sprintf('\n%.1f%%', res$convProbBbeatsA*100),
          sprintf('\n%.1f%%', res$convProbCbeatsA*100),
          sprintf('\n%.1f%%', res$convProbCbeatsB*100)
        ),
        ARPPU = c(
          sprintf('\n%.1f%%', res$revProbBbeatsA*100),
          sprintf('\n%.1f%%', res$revProbCbeatsA*100),
          sprintf('\n%.1f%%', res$revProbCbeatsB*100)
        ),
        ARPU = c(
          sprintf('\n%.1f%%', res$arpuProbBbeatsA*100),
          sprintf('\n%.1f%%', res$arpuProbCbeatsA*100),
          sprintf('\n%.1f%%', res$arpuProbCbeatsB*100)
        ),
        Retention = c(
          sprintf('\n%.2g%% [\n%.2g%%]',as.numeric((convret_B-convret_A)/convret_A)*100, prob_B_beats_A(alpharet_A, betaret_A, alpharet_B, betaret_B)*100),
          sprintf('\n%.2g%% [\n%.2g%%]',as.numeric((convret_C-convret_A)/convret_A)*100, prob_C_beats_A(alpharet_A, betaret_A, alpharet_C, betaret_C)*100),
          sprintf('\n%.2g%% [\n%.2g%%]',as.numeric((convret_C-convret_B)/convret_B)*100, prob_C_beats_B(alpharet_B, betaret_B, alpharet_C, betaret_C)*100)
          
        )
      )
      colnames(tab) <- c(' ', 'Conversion', 'ARPPU', 'ARPU', 'Retention')
      tab
    }, spacing = 'xs')
    
    output$table3 <- renderTable({
      tab <- data.frame(
        column1 = c(
          'Probability that A is better than B',
          'Probability that A is better than C',
          'Probability that B is better than C'
        ),
        conversion = c(
          sprintf('\n%.1f%%', res$convProbAbeatsB*100),
          sprintf('\n%.1f%%', res$convProbAbeatsC*100),
          sprintf('\n%.1f%%', res$convProbBbeatsC*100)
        ),
        ARPPU = c(
          sprintf('\n%.1f%%', res$revProbAbeatsB*100),
          sprintf('\n%.1f%%', res$revProbAbeatsC*100),
          sprintf('\n%.1f%%', res$revProbBbeatsC*100)
        ),
        ARPU = c(
          sprintf('\n%.1f%%', res$arpuProbAbeatsB*100),
          sprintf('\n%.1f%%', res$arpuProbAbeatsC*100),
          sprintf('\n%.1f%%', res$arpuProbBbeatsC*100)
        ),
        Retention1 = c(
          sprintf('\n%.2g%% [\n%.2g%%]', as.numeric((convret_A-convret_B)/convret_B)*100, prob_A_beats_B(alpharet_B, betaret_B, alpharet_A, betaret_A)*100),
          sprintf('\n%.2g%% [\n%.2g%%]', as.numeric((convret_A-convret_C)/convret_C)*100, prob_A_beats_C(alpharet_C, betaret_C, alpharet_A, betaret_A)*100),
          sprintf('\n%.2g%% [\n%.2g%%]', as.numeric((convret_B-convret_C)/convret_C)*100, prob_B_beats_C(alpharet_C, betaret_C, alpharet_B, betaret_B)*100)
        )
      )
      colnames(tab) <- c(' ', 'Conversion', 'ARPPU', 'ARPU', 'Retention')
      tab
    }, spacing = 'xs')
    
    output$table4 <- renderTable({
      tab <- data.frame(
        column1 = c(
          '*Retention1 refers to A better than B, A better than C, B better than C',
          '*Retention2 refers to B better than A, C better than A, C better than B',
          'Report largest value (positive integer)'
        )
      )
      colnames(tab) <- c(' ')
      tab
    }, spacing = 'xs')
    
    output$table5 <- renderTable({
      tab <- data.frame(
        better = c(
          paste('On ', input$testday, ', the Control group had a ', round(res$convProbAbeatsB*100, digits=2), '% greater conversion rate than Group B, and a ', round(res$convProbAbeatsC*100, digits=2),'% greater conversion rate compared to Group C.', sep=''),
          paste('On ', input$testday, ', Group B had a ', round(res$convProbBbeatsA*100, digits=2), '% greater conversion rate compared to the Control group, and a ', round(res$convProbBbeatsC*100, digits=2),'% greater conversion rate compared to Group C.', sep=''),
          paste('On ', input$testday, ', Group C had a ', round(res$convProbCbeatsA*100, digits=2), '% greater conversion rate compared to the Control group, and a ', round(res$convProbCbeatsB*100, digits=2),'% greater conversion rate compared to Group B.', sep='')
          
        )
      )
      colnames(tab) <- c(' ')
      tab
    }, spacing = 'xs')
    
    
    output$table6 <- renderTable({
      tab <- data.frame(
        better = c(
          paste('On ', input$testday, ', the Control group had a ', round(res$arpuProbAbeatsB*100,digits=2), '% greater lifetime value than Group B, and a ', round(res$arpuProbAbeatsC*100,digits=2), '% greater lifetime value compared to Group C.', sep=''),
          paste('On ', input$testday, ', Group B had a ', round(res$arpuProbBbeatsA*100, digits=2), '% greater lifetime value compared to the Control group, and a ', round(res$arpuProbBbeatsC*100, digits=2), '% greater lifetime value compared to Group C.', sep=''),
          paste('On ', input$testday, ', Group C had a ', round(res$arpuProbCbeatsA*100, digits=2), '% greater lifetime value compared to the Control group, and a ', round(res$arpuProbCbeatsB*100, digits=2), '% greater lifetime value compared to Group B.', sep='')
        )
      )
      colnames(tab) <- c(' ')
      tab
    }, spacing = 'xs')
    
    
    output$table7 <- renderTable({
      tab <- data.frame(
        better = c(
          paste('On ', input$testday,', the Control group had a ', round(((convret_A-convret_B)/convret_B)*100,digits=2), '% greater retention rate compared to Group B (', round(prob_A_beats_B(alpharet_B, betaret_B, alpharet_A, betaret_A)*100,digits =0),
                '% confidence), and a ', round(((convret_A-convret_C)/convret_C)*100,digits=2),'% greater retention rate compared to Group C (', round(prob_A_beats_C(alpharet_C, betaret_C, alpharet_A, betaret_A)*100,digits=0),
                '% confidence).', sep=''),
          
          paste('On ', input$testday,', Group B had a ', round(((convret_B-convret_A)/convret_A)*100,digits=2), '% greater retention rate compared to the Control group (', round(prob_B_beats_A(alpharet_A, betaret_A, alpharet_B, betaret_B)*100,digits =0),
                '% confidence), and a ', round(((convret_B-convret_C)/convret_C)*100,digits=2),'% greater retention rate compared to Group C (', round(prob_B_beats_C(alpharet_C, betaret_C, alpharet_B, betaret_B)*100,digits =0),
                '% confidence).', sep=''),
          
          paste('On ', input$testday,', Group C had a ', round(((convret_C-convret_A)/convret_A)*100,digits=2), '% greater retention rate compared to the Control group (', round(prob_C_beats_A(alpharet_A, betaret_A, alpharet_C, betaret_C)*100,digits =0),
                '% confidence), and a ', round(((convret_C-convret_B)/convret_B)*100,digits=2),'% greater retention rate compared to Group B (', round(prob_C_beats_B(alpharet_B, betaret_B, alpharet_C, betaret_C)*100,digits =0),
                '% confidence).', sep='')
        )
      )
      colnames(tab) <- c(' ')
      tab
    }, spacing = 'xs')
  })
})
